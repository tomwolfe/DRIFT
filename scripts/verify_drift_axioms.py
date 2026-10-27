#!/usr/bin/env python3
"""Certify the axiom set of the four exported DRIFT Hill lemmas.

``QED/verify_pbpk_lemmas.py`` proves the lemmas and rejects ``sorry``, but it
does not expose the ``#print axioms`` set.  This script closes that gap: it
emits the four theorems with the same tactics QED found, compiles the file
with the pinned Lean toolchain, parses every ``#print axioms`` line and
requires the set to be exactly ``[propext, Classical.choice, Quot.sound]``.

Writes ``/tmp/drift_axioms.json`` for the tether evidence collector.

Usage:
    python3 scripts/verify_drift_axioms.py [--qed-dir ../QED]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

EXPECTED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}

HEADER = """import Mathlib.Tactic
import Mathlib.Data.Real.Basic
"""


def _theorem(name: str, statement: str, tactic: str) -> str:
    return f"theorem {name} {statement} := by\n  {tactic}\n"


# These are the declarations QED proved for /tmp/drift_lemmas.txt.  Keeping the
# proofs here (rather than re-deriving them) is deliberate: the point of this
# gate is the AXIOM SET, not a second tactic search.
THEOREMS = [
    _theorem(
        "hill_nonneg",
        "(C : ℝ) (Kd : ℝ) (n : ℝ) (hC : 0 < C) (hKd : 0 < Kd) (hn : 0 < n) :\n"
        "    C ^ n / (C ^ n + Kd ^ n) >= 0",
        "intros; positivity",
    ),
    _theorem(
        "hill_le_one",
        "(C : ℝ) (Kd : ℝ) (n : ℝ) (hC : 0 < C) (hKd : 0 < Kd) (hn : 0 < n) :\n"
        "    C ^ n / (C ^ n + Kd ^ n) + Kd ^ n / (C ^ n + Kd ^ n) = 1",
        "field_simp",
    ),
    _theorem(
        "hill_mono",
        "(C1 : ℝ) (C2 : ℝ) (Kd : ℝ) (n : ℝ) (hC1 : 0 < C1) (hC2 : 0 < C2)\n"
        "    (hKd : 0 < Kd) (hn : 0 < n) :\n"
        "    (C2 ^ n / (C2 ^ n + Kd ^ n) - C1 ^ n / (C1 ^ n + Kd ^ n))"
        " * (C2 ^ n + Kd ^ n) * (C1 ^ n + Kd ^ n)"
        " = (C2 ^ n - C1 ^ n) * Kd ^ n",
        "intros; field_simp; ring",
    ),
    _theorem(
        "hill_limits",
        "(C : ℝ) (Kd : ℝ) (n : ℝ) (hC : 0 < C) (hKd : 0 < Kd) (hn : 0 < n) :\n"
        "    C ^ n / (C ^ n + Kd ^ n) = 1 / (1 + Kd ^ n / C ^ n)",
        "field_simp",
    ),
]


def _lean_path(qed: Path) -> str:
    parts = [str(qed / ".lake" / "build" / "lib" / "lean")]
    packages = qed / ".lake" / "packages"
    if packages.is_dir():
        for pkg in sorted(packages.iterdir()):
            cand = pkg / ".lake" / "build" / "lib" / "lean"
            if cand.is_dir():
                parts.append(str(cand))
    return os.pathsep.join(parts)


def _lean_bin() -> list[str]:
    found = shutil.which("lean")
    if found:
        return [found]
    elan = shutil.which("elan")
    if elan:
        return [elan, "show", "lean"]
    return ["lean"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--qed-dir", type=Path,
                    default=Path(__file__).resolve().parents[2] / "QED")
    args = ap.parse_args(argv)
    qed = args.qed_dir.resolve()
    if not (qed / "lakefile.lean").is_file():
        print(f"QED checkout not found at {qed}", file=sys.stderr)
        return 1

    body = HEADER + "\n".join(THEOREMS)
    body += "\n" + "\n".join(f"#print axioms {name}" for name in
                             ("hill_nonneg", "hill_le_one", "hill_mono",
                              "hill_limits")) + "\n"
    tmp = qed / "_drift_axiom_check.lean"
    tmp.write_text(body, encoding="utf-8")
    env = dict(os.environ)
    env["LEAN_PATH"] = _lean_path(qed)
    try:
        proc = subprocess.run(
            [*_lean_bin(), str(tmp)], cwd=str(qed), capture_output=True,
            text=True, timeout=600, env=env)
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"DRIFT axiom check failed to run: {exc}", file=sys.stderr)
        return 1
    finally:
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass

    output = proc.stdout + proc.stderr
    if proc.returncode != 0:
        print("DRIFT axiom check: Lean compilation failed", file=sys.stderr)
        print(output[-2000:], file=sys.stderr)
        return 1
    if "sorry" in output or "sorryAx" in output:
        print("DRIFT axiom check: sorry lemma in output", file=sys.stderr)
        return 1
    groups = re.findall(r"depends on axioms:?\s*\[(.*?)\]", output)
    printed = len(re.findall(r"depends on axioms", output))
    if not groups or len(groups) < printed:
        print(
            f"DRIFT axiom check: could not parse #print axioms "
            f"({len(groups)} of {printed})", file=sys.stderr)
        return 1
    axioms = sorted({a.strip().strip("'") for grp in groups
                     for a in grp.split(",") if a.strip()})
    if set(axioms) != EXPECTED_AXIOMS:
        print(
            f"DRIFT axiom check: unexpected axioms {axioms}; expected "
            f"{sorted(EXPECTED_AXIOMS)}", file=sys.stderr)
        return 1

    payload = {
        "lemmas": ["hill_nonneg", "hill_le_one", "hill_mono", "hill_limits"],
        "axioms": axioms,
        "verified": True,
    }
    Path("/tmp/drift_axioms.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"DRIFT axioms certified: {axioms}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
