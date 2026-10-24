#!/usr/bin/env python3
"""Verify DRIFT's exported PD/metabolic lemmas through QED.

This is the DRIFT analogue of ``VeriTrial/scripts/verify_formal_gate.py``, and
it is deliberately placed on the DRIFT side: QED is domain-agnostic, and a
file named after a domain inside ``QED/scripts/`` would be exactly the
leakage the QED repository refuses.  The VERIFIER is QED's own, generic one
(``QED/verify_pbpk_lemmas.py``), driven from here.

The gate is fail-closed on four counts, and each of them is a way a green
result could be hollow:

  1. The supplied lemma file must be EXACTLY what DRIFT's own exporter
     emits right now.  A hand-maintained or stale file is refused, so the
     gate cannot certify a set the model does not produce.
  2. No lemma may be THEATER: a reflexive identity (``X = X``) or a
     closed numeric statement (``129 = 129``) proves nothing about the
     model, so such a lemma is refused rather than certified.
  3. No lemma may contain a ``sorry``/``sorryAx`` placeholder.
  4. QED must verify every lemma; one failure fails the gate.

THE CERTIFICATION IS PARALLEL, AND THE PARTITION IS CHECKED
-----------------------------------------------------------
Proving one lemma takes QED several seconds (each tactic attempt is a Lean
compile), so proving the set one lemma at a time makes the gate the most
expensive command in the stack -- and the mutation meta-check re-runs this
gate once per mutant, which does not fit the mission's wall budget.  The
certification is therefore split into DISJOINT CHUNKS that QED's own generic
verifier proves CONCURRENTLY.

The split changes the SCHEDULING, not the certificate:

  * the chunks are a partition of the lemma set -- every lemma is proved
    once and none twice, and the partition is asserted before a single
    proof is attempted;
  * every chunk must exit zero, i.e. QED's own verifier must accept every
    lemma in it;
  * the ``verified:`` lines of all chunks together must account for the
    lemma set EXACTLY -- a lemma that is missing, or proved twice, or
    silently dropped fails the gate;
  * a skipped (Mathlib-dependent) lemma still fails the gate under
    ``--strict``.

So the verdict is the same single verdict the sequential run produced, and
it is checked more strictly: the sequential run only asked QED's exit code,
this one also asks which lemmas were proved.

Usage:
    python3 scripts/verify_drift_lemmas.py <lemmas_file> [--strict]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

DRIFT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = DRIFT_ROOT / "scripts"

#: How many lemmas one QED process proves.  The chunk size is the knob that
#: trades concurrency against memory: each lemma is a Lean compile, so a
#: chunk of four keeps the machine busy without asking Lean for more
#: processes than it has cores.
CHUNK_SIZE = 4


def qed_dir() -> Path:
    """Locate the QED checkout (sibling of this repository)."""
    env = sys.modules["os"].environ.get("QED_DIR")
    if env:
        return Path(env).resolve()
    return (DRIFT_ROOT.parent / "QED").resolve()


def live_lemmas() -> list[str]:
    """The lemmas DRIFT's exporter emits right now."""
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    import export_drift_to_qed as ex  # noqa: PLC0415
    return [l for l in ex.build_lemmas(parametric=True)
            if l.strip() and not l.strip().startswith("--")]


def _is_sorry_placeholder(lemma: str) -> bool:
    return bool(re.search(r"\bsorry\b|\bsorryAx\b", lemma))


def _split_sum(text: str, operator: str) -> list[str]:
    """Split on a top-level operator, respecting parentheses."""
    parts: list[str] = []
    depth = 0
    start = 0
    for index, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == operator and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def _is_commutative_tautology(left: str, right: str) -> bool:
    """True when the two sides are the same terms in a different order.

    ``9 * p + 1 = 1 + 9 * p`` and ``a * b = b * a`` are true for every
    model, so they certify nothing.  A purely textual check misses them;
    comparing the top-level summands (or factors) as multisets catches the
    reordering without needing a symbolic engine.
    """
    for operator in ("+", "*"):
        lhs = _split_sum(left, operator)
        rhs = _split_sum(right, operator)
        if len(lhs) > 1 and len(rhs) > 1 and sorted(lhs) == sorted(rhs):
            return True
    return False


def _is_theater(lemma: str) -> bool:
    """True for a lemma that proves nothing about the model.

    Three classes, all refused:

    * a reflexive identity (``X = X``) -- the statement is true whatever
      the model does;
    * a closed numeric statement (``129 = 129``) -- no variable at all;
    * a commutative tautology (``9 * p + 1 = 1 + 9 * p``) -- the same
      terms in a different order, true for every model.
    """
    text = lemma.strip()
    if "=" not in text:
        return True
    left, _, right = text.partition("=")
    if left.strip() == right.strip():
        return True
    if not re.search(r"[A-Za-z_]", left):
        return True
    if _is_commutative_tautology(left.strip(), right.strip()):
        return True
    return False


def chunk_lemmas(lemmas: list[str],
                 chunk_size: int = CHUNK_SIZE) -> list[list[str]]:
    """Split the lemma set into DISJOINT chunks, in order.

    The chunks are a partition: every lemma is in exactly one chunk, so the
    certification proves the whole set once and nothing twice.  The caller
    asserts the partition before a single proof is attempted.
    """
    size = max(1, chunk_size)
    return [lemmas[i:i + size] for i in range(0, len(lemmas), size)]


def _verified_lemmas(stdout: str) -> list[str]:
    """The lemmas QED's verifier reported as proved."""
    return [line[len("verified: "):].strip()
            for line in stdout.splitlines()
            if line.startswith("verified: ")]


def _skipped_lemmas(stdout: str) -> list[str]:
    """The lemmas QED's verifier skipped (a Mathlib-dependent lemma)."""
    return [line.strip() for line in stdout.splitlines()
            if line.startswith("skipped")]


@contextmanager
def _preserved_qed_artifact(qed: Path):
    """Put QED's tracked artefact back byte for byte after the run.

    QED's own verifier already snapshots and restores ``output.lean``, but
    the chunk processes run CONCURRENTLY, so the artefact's final bytes would
    otherwise depend on how their writes interleaved.  The gate must not
    leave a mark on the prover's repository, so the whole parallel phase is
    wrapped in a snapshot/restore of the same bytes: the prover's tree comes
    out exactly as it went in, whatever the interleaving was.
    """
    artifact = qed / "output.lean"
    try:
        original = artifact.read_bytes()
        existed = True
    except OSError:
        original = None
        existed = False
    try:
        yield
    finally:
        try:
            if existed and original is not None:
                artifact.write_bytes(original)
            elif not existed:
                artifact.unlink(missing_ok=True)
        except OSError:
            pass


def _prove_chunks(verifier: Path, qed: Path,
                  chunk_files: list[Path]) -> list[subprocess.CompletedProcess]:
    """Run QED's generic verifier on every chunk, concurrently."""
    def _prove(chunk_file: Path) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(verifier), str(chunk_file)],
            cwd=str(qed), capture_output=True, text=True, timeout=7200,
        )

    if len(chunk_files) == 1:
        return [_prove(chunk_files[0])]
    with ThreadPoolExecutor(max_workers=len(chunk_files)) as pool:
        return list(pool.map(_prove, chunk_files))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("lemmas", type=Path,
                    help="file of DRIFT lemmas, one per line")
    ap.add_argument("--strict", action="store_true",
                    help="fail if a Mathlib-dependent lemma is skipped")
    args = ap.parse_args(argv)

    if not args.lemmas.is_file():
        print(f"DRIFT FORMAL GATE FAILED: lemmas file not found: "
              f"{args.lemmas}", file=sys.stderr)
        return 1

    file_lemmas = [
        line.strip()
        for line in args.lemmas.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("--")
    ]
    if not file_lemmas:
        print("DRIFT FORMAL GATE FAILED: the lemma file is empty",
              file=sys.stderr)
        return 1

    # 1. Single source of truth: the file must be what DRIFT emits now.
    emitted = live_lemmas()
    if set(file_lemmas) != set(emitted):
        missing = sorted(set(emitted) - set(file_lemmas))
        extra = sorted(set(file_lemmas) - set(emitted))
        print("DRIFT FORMAL GATE FAILED (fail-closed): the lemma file is "
              "not the set the DRIFT source emits.", file=sys.stderr)
        if missing:
            print(f"  missing from file: {missing}", file=sys.stderr)
        if extra:
            print(f"  not produced by the DRIFT exporter: {extra}",
                  file=sys.stderr)
        return 1

    # 2. No theater.
    theater = [l for l in file_lemmas if _is_theater(l)]
    if theater:
        print("DRIFT FORMAL GATE FAILED (fail-closed): the lemma set "
              "contains statements that prove nothing about the model:",
              file=sys.stderr)
        for lemma in theater:
            print(f"  {lemma}", file=sys.stderr)
        return 1

    # 3. No sorry, before QED is even asked.
    sorry = [l for l in file_lemmas if _is_sorry_placeholder(l)]
    if sorry:
        print("DRIFT FORMAL GATE FAILED: sorry placeholder in the lemma set",
              file=sys.stderr)
        return 1

    qed = qed_dir()
    verifier = qed / "verify_pbpk_lemmas.py"
    if not qed.is_dir() or not verifier.is_file():
        print(f"DRIFT FORMAL GATE FAILED (fail-closed): QED not found at "
              f"{qed}; set QED_DIR or place QED as a sibling of this "
              "repository.", file=sys.stderr)
        return 1

    # The certification is split into chunks that are proved concurrently.
    # The split must be a PARTITION: nothing dropped, nothing proved twice.
    chunks = chunk_lemmas(file_lemmas)
    flattened = [lemma for chunk in chunks for lemma in chunk]
    if flattened != file_lemmas or len(set(flattened)) != len(flattened):
        print("DRIFT FORMAL GATE FAILED (fail-closed): the certification "
              "chunks are not a partition of the lemma set.", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory(prefix="drift-formal-gate-") as tmp:
        chunk_files: list[Path] = []
        for index, chunk in enumerate(chunks):
            path = Path(tmp) / f"chunk-{index:02d}.txt"
            path.write_text("\n".join(chunk) + "\n", encoding="utf-8")
            chunk_files.append(path)

        with _preserved_qed_artifact(qed):
            results = _prove_chunks(verifier, qed, chunk_files)

    verified: list[str] = []
    skipped: list[str] = []
    failed = False
    for proc in results:
        sys.stdout.write(proc.stdout)
        if proc.stderr:
            sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            failed = True
        verified.extend(_verified_lemmas(proc.stdout))
        skipped.extend(_skipped_lemmas(proc.stdout))

    if failed:
        print("DRIFT FORMAL GATE FAILED: QED did not verify every lemma.",
              file=sys.stderr)
        return 1

    # Every lemma must be accounted for EXACTLY once: a lemma that is
    # missing from the verdicts, or counted twice, fails the gate.
    if sorted(verified) != sorted(file_lemmas) \
            or len(verified) != len(set(verified)):
        missing = sorted(set(file_lemmas) - set(verified))
        extra = sorted(set(verified) - set(file_lemmas))
        print("DRIFT FORMAL GATE FAILED (fail-closed): the certification did "
              "not account for every lemma exactly once.", file=sys.stderr)
        if missing:
            print(f"  not verified: {missing}", file=sys.stderr)
        if extra:
            print(f"  verified but not in the lemma set: {extra}",
                  file=sys.stderr)
        return 1

    if skipped and args.strict:
        print("DRIFT FORMAL GATE FAILED (--strict): a Mathlib-dependent "
              "lemma was skipped, so the gate would degrade silently.",
              file=sys.stderr)
        for line in skipped:
            print(f"  {line}", file=sys.stderr)
        return 1

    print(f"DRIFT FORMAL GATE PASSED: all {len(file_lemmas)} lemmas "
          "verified by QED (no sorry).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
