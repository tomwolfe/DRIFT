"""The formal export and the formal gate's own fail-closed checks.

The export and the verifier are the gate's own code, so they need their own
teeth: a mutation of the export (say, one that stops reading the source)
would otherwise certify a model that no longer exists, and a mutation of
the verifier's THEATER filter would let a reflexive identity through.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

DRIFT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = DRIFT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import export_drift_to_qed as ex  # noqa: E402
import verify_drift_lemmas as vdl  # noqa: E402


def test_the_export_reads_the_source_and_fails_closed(monkeypatch) -> None:
    """A model that no longer implements the Hill equation is REFUSED."""
    real = ex._binding_source

    def _broken() -> Path:
        path = DRIFT_ROOT / "tests" / "_tmp_broken_binding.py"
        path.write_text("class BindingEngine:\n    pass\n", encoding="utf-8")
        return path

    monkeypatch.setattr(ex, "_binding_source", _broken)
    try:
        with pytest.raises(ValueError):
            ex.build_hill_lemmas(ex.extract_hill_shape())
    finally:
        (DRIFT_ROOT / "tests" / "_tmp_broken_binding.py").unlink(
            missing_ok=True)
        monkeypatch.setattr(ex, "_binding_source", real)


def test_the_export_refuses_a_changed_flux_or_bridge(monkeypatch) -> None:
    """The other two certified equations are checked the same way."""
    for target, reader in (("flux", ex.extract_flux_shape),
                            ("bridge", ex.extract_bridge_scaling),
                            ("milstein", ex.extract_milstein_step)):
        real = ex._metabolic_source if target in ("flux", "bridge") \
            else ex._signaling_source
        path = DRIFT_ROOT / "tests" / f"_tmp_broken_{target}.py"
        path.write_text("# nothing here\n", encoding="utf-8")
        monkeypatch.setattr(ex, "_metabolic_source", lambda p=path: p)
        monkeypatch.setattr(ex, "_signaling_source", lambda p=path: p)
        try:
            with pytest.raises(ValueError):
                reader()
        finally:
            path.unlink(missing_ok=True)
            monkeypatch.setattr(ex, "_metabolic_source", real)
            monkeypatch.setattr(ex, "_signaling_source", real)


def test_the_export_refuses_a_changed_inhibition_modulation(monkeypatch) -> None:
    """``drift/topology.py``'s inhibition algebra is pinned, not asserted.

    Without this reader the three inhibition lemmas would be hand-written
    statements the source is never checked against -- a changed modulation
    could keep the certificate.
    """
    path = DRIFT_ROOT / "tests" / "_tmp_changed_topology.py"
    path.write_text(
        "def pi3k_akt_mtor_drift(state, params, feedback=None):\n"
        "    effective_pi3k = pi3k / (1.0 - inh_pi3k)\n"
        "    effective_akt = akt * (1.0 - inh_akt)\n"
        "    dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - 1\n"
        "    return dakt\n", encoding="utf-8")
    real = ex._topology_source
    monkeypatch.setattr(ex, "_topology_source", lambda: path)
    try:
        with pytest.raises(ValueError):
            ex.extract_inhibition_shape()
    finally:
        path.unlink(missing_ok=True)
        monkeypatch.setattr(ex, "_topology_source", real)


def test_the_export_reads_the_ast_not_the_text(tmp_path) -> None:
    """Reformatting is the same equation; a changed operator is not."""
    source = (DRIFT_ROOT / "drift" / "binding.py").read_text(
        encoding="utf-8")

    reformatted = tmp_path / "reformatted.py"
    reformatted.write_text(
        source.replace(
            "numerator = drug_concentration**self.hill_coefficient",
            "numerator = (drug_concentration  **  self.hill_coefficient)"),
        encoding="utf-8")
    assert ex.extract_hill_shape(reformatted)["concentration"] == "C_drug", \
        "a reformatted equation was refused, so the reader is textual"

    changed = tmp_path / "changed.py"
    changed.write_text(
        source.replace(
            "denominator = numerator + kd**self.hill_coefficient",
            "denominator = numerator + kd*self.hill_coefficient"),
        encoding="utf-8")
    with pytest.raises(ValueError):
        ex.extract_hill_shape(changed)


def test_the_ast_reader_selects_the_requested_annotated_assignment() -> None:
    """``_assigned`` must read the name it was asked for.

    The annotated-assignment branch has to compare the target's name;
    accepting ANY annotated target would return the first, unrelated
    value in the function and certify an expression the model does not
    contain.  Two annotated assignments make the difference visible.
    """
    func = ast.parse(
        "def f():\n"
        "    unrelated: int = 1\n"
        "    wanted: int = 2\n"
        "    return wanted\n"
    ).body[0]

    node = ex._assigned(func, "wanted", "tests")

    assert ast.dump(node) == ast.dump(ast.parse("2", mode="eval").body)


def test_the_signaling_builder_refuses_a_missing_inhibition_shape() -> None:
    """The inhibition lemmas are tied to the shape the reader returned.

    ``extract_inhibition_shape`` verifies the three source expressions and
    returns the modulation it read.  A reader that returns nothing (or a
    different modulation) must fail the export closed instead of emitting
    lemmas with no source behind them.
    """
    real = ex.extract_inhibition_shape()
    assert ex.build_signaling_lemmas(real)  # the real shape is accepted
    with pytest.raises(ValueError):
        ex.build_signaling_lemmas(None)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        ex.build_signaling_lemmas(
            {"species": "akt", "inhibition": "inh_akt"})


def test_the_exported_lemmas_have_no_decimal_literal() -> None:
    """A decimal literal does not tokenize, so it would fail closed."""
    import re
    for lemma in ex.build_lemmas(parametric=True):
        assert not re.search(r"\d\.\d", lemma), lemma


def test_the_theater_filter_refuses_a_bare_expression() -> None:
    """A statement with no '=' is not a lemma, so it is refused."""
    assert vdl._is_theater("x + 1") is True
    assert vdl._is_theater("sorry") is True


def test_the_gate_preserves_the_qed_artifact(tmp_path) -> None:
    """The parallel certification must leave QED's output.lean byte for
    byte, whether it existed before the run or not."""
    existing = tmp_path / "with_artifact"
    existing.mkdir()
    artifact = existing / "output.lean"
    artifact.write_bytes(b"-- original certificate\n")
    with vdl._preserved_qed_artifact(existing):
        artifact.write_bytes(b"-- clobbered by a chunk\n")
    assert artifact.read_bytes() == b"-- original certificate\n", \
        "the gate did not restore QED's pre-existing output.lean"

    absent = tmp_path / "without_artifact"
    absent.mkdir()
    with vdl._preserved_qed_artifact(absent):
        (absent / "output.lean").write_bytes(b"-- created by a chunk\n")
    assert not (absent / "output.lean").exists(), \
        "the gate left a new output.lean behind"


def test_the_theater_filter_refuses_a_reflexive_and_a_numeric_lemma() -> None:
    """A lemma that proves nothing about the model is refused."""
    assert vdl._is_theater("X = X") is True
    assert vdl._is_theater("129 = 129") is True
    assert vdl._is_theater("0 = 0") is True
    assert vdl._is_theater(
        "9 * protein_level + 1 = 1 + 9 * protein_level") is True
    assert vdl._is_theater("a * b = b * a") is True
    assert vdl._is_theater("C_drug / (C_drug + Kd) >= 0") is False
    assert vdl._is_theater(
        "C_drug / (C_drug + Kd) + Kd / (C_drug + Kd) = 1") is False


def test_the_theater_filter_is_not_vacuous() -> None:
    """A filter that refuses everything is not a filter."""
    real = ex.build_lemmas(parametric=True)
    assert real, "the export emitted no lemmas"
    assert not any(vdl._is_theater(lemma) for lemma in real), (
        "the export emitted a lemma the gate would refuse: the two "
        "components disagree")


def test_the_spec_blocks_match_the_exported_lemmas() -> None:
    """Each .tex spec carries exactly the lemmas the exporter emits."""
    blocks = ex.build_spec_blocks()
    assert blocks, "no spec blocks were built"
    for stem, lemmas in blocks.items():
        path = ex.spec_dir() / f"{stem}.tex"
        assert path.is_file(), path
        text = path.read_text(encoding="utf-8")
        assert ex.SPEC_BLOCK_BEGIN in text and ex.SPEC_BLOCK_END in text
        body = text.split(ex.SPEC_BLOCK_BEGIN, 1)[1].split(
            ex.SPEC_BLOCK_END, 1)[0]
        stated = [line[len("%% QED "):].strip()
                  for line in body.splitlines()
                  if line.strip().startswith("%% QED ")]
        assert stated == lemmas, (stem, stated, lemmas)


def test_the_gate_refuses_a_lemma_file_that_is_not_the_export() -> None:
    """A hand-edited or stale lemma file is refused, not certified."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "lemmas.txt"
        path.write_text("C_drug / (C_drug + Kd) >= 0\n", encoding="utf-8")
        rc = vdl.main([str(path)])
        assert rc == 1, "the gate accepted a lemma file that is not the export"


def test_the_export_cli_fails_closed_on_a_changed_model(
        monkeypatch, tmp_path) -> None:
    """The CLI must EXIT NON-ZERO when the source stops implementing the
    certified equation -- a ValueError swallowed into a zero exit would
    look like a green export."""
    def _broken(parametric: bool = True):
        raise ValueError("source no longer matches the certified equation")

    monkeypatch.setattr(ex, "build_lemmas", _broken)
    out = tmp_path / "lemmas.txt"
    assert ex.main(["--out", str(out)]) == 1, \
        "the export CLI returned success for a fail-closed error"


def test_the_flux_lemmas_state_the_implemented_normalized_mapping() -> None:
    """The code returns C*(Km+1)/(Km+C), so the certificate must too.

    Certifying the textbook C/(Km+C) would leave the implemented mapping
    uncovered even though the prover is happy -- the review gate caught
    exactly that.  This test pins the factor and rejects the old unit
    budget lemma, which is false for the normalized mapping at large C.
    """
    lemmas = ex.build_flux_lemmas(ex.extract_flux_shape())
    assert lemmas, "the flux spec emitted no lemmas"
    assert any("(Km + 1)" in lemma for lemma in lemmas), \
        "the implemented (Km + 1) normalization is missing from the lemmas"
    assert not any(lemma.startswith("C_flux / (Km + C_flux) + Km /")
                   for lemma in lemmas), \
        "the textbook unit-budget lemma survived; it is not this mapping"


def test_the_hill_lemmas_carry_the_source_exponent(tmp_path) -> None:
    """The certificate is specialised to the default hill_coefficient, so
    changing that default changes the lemma text instead of silently
    keeping yesterday's certificate."""
    source = (DRIFT_ROOT / "drift" / "binding.py").read_text(encoding="utf-8")
    two = tmp_path / "binding_two.py"
    two.write_text(source.replace(
        "hill_coefficient: float = 1.0", "hill_coefficient: float = 2.0"),
        encoding="utf-8")
    shape = ex.extract_hill_shape(two)
    assert shape["exponent_value"] == "2"
    lemmas = ex.build_hill_lemmas(shape)
    assert all("^2" in lemma for lemma in lemmas if "C_drug" in lemma), \
        f"the exponent was not carried into the lemmas: {lemmas}"
    assert lemmas != ex.build_hill_lemmas(ex.extract_hill_shape()), \
        "a changed exponent produced the same certificate"

    fractional = tmp_path / "binding_fractional.py"
    fractional.write_text(source.replace(
        "hill_coefficient: float = 1.0", "hill_coefficient: float = 1.5"),
        encoding="utf-8")
    with pytest.raises(ValueError):
        ex.extract_hill_shape(fractional)


def test_the_gate_refuses_a_missing_lemma_file(tmp_path) -> None:
    """A missing lemma file is a failed gate, not a vacuous pass."""
    missing = tmp_path / "does_not_exist.txt"
    assert vdl.main([str(missing)]) == 1, \
        "the gate did not fail when the lemma file was missing"


def test_the_gate_refuses_a_non_partition_of_the_lemma_set(
        monkeypatch, tmp_path) -> None:
    """A chunking bug that drops a lemma fails the gate before proving."""
    lemma = "C_drug / (C_drug + Kd) >= 0"
    path = tmp_path / "lemmas.txt"
    path.write_text(lemma + "\n", encoding="utf-8")
    monkeypatch.setattr(vdl, "live_lemmas", lambda: [lemma])
    monkeypatch.setattr(vdl, "chunk_lemmas", lambda lemmas: [[]])
    assert vdl.main([str(path)]) == 1, \
        "the gate certified a chunking that dropped every lemma"


def test_the_gate_refuses_a_strict_skip(monkeypatch, tmp_path) -> None:
    """Under --strict a skipped lemma fails the gate, even if the rest
    verified; otherwise the certification would degrade silently."""
    import subprocess

    lemma = "C_drug / (C_drug + Kd) >= 0"
    path = tmp_path / "lemmas.txt"
    path.write_text(lemma + "\n", encoding="utf-8")
    monkeypatch.setattr(vdl, "live_lemmas", lambda: [lemma])

    def _fake_prove(*_args, **_kwargs):
        return [subprocess.CompletedProcess(
            args=[], returncode=0,
            stdout=f"verified: {lemma}\nskipped: {lemma} (needs Mathlib)\n",
            stderr="")]

    monkeypatch.setattr(vdl, "_prove_chunks", _fake_prove)
    assert vdl.main([str(path), "--strict"]) == 1, \
        "a skipped lemma was certified under --strict"
    assert vdl.main([str(path)]) == 0, \
        "a skipped lemma was allowed to fail a non-strict gate"


def test_the_certification_chunks_partition_the_lemma_set() -> None:
    """No lemma is dropped or proved twice when the gate chunks the set."""
    lemmas = ex.build_lemmas(parametric=True)
    chunks = vdl.chunk_lemmas(lemmas)
    assert chunks, "the lemma set produced no chunk"
    assert all(chunk for chunk in chunks), "an empty chunk would waste a prover"
    flat = [lemma for chunk in chunks for lemma in chunk]
    assert flat == lemmas, "the chunks are not the lemma set, in order"
    assert len(set(flat)) == len(flat), "a lemma is in more than one chunk"


def test_the_spec_sync_refuses_a_file_without_a_machine_checkable_block(
        tmp_path) -> None:
    """A spec with no machine-checkable block is refused, not rewritten."""
    path = tmp_path / "no_block.tex"
    path.write_text("% prose only\n", encoding="utf-8")
    with pytest.raises(ValueError):
        ex.sync_spec_file(path, ["X = 1"])


def test_the_spec_sync_says_whether_it_rewrote_the_block(tmp_path) -> None:
    """The sync reports True/False -- a stale block is rewritten, a fresh
    one is left alone."""
    lemmas = ex.build_hill_lemmas(ex.extract_hill_shape())
    path = tmp_path / "spec.tex"
    path.write_text(
        f"% prose\n{ex.SPEC_BLOCK_BEGIN}\n% stale\n{ex.SPEC_BLOCK_END}\n",
        encoding="utf-8")
    assert ex.sync_spec_file(path, lemmas) is True, \
        "the stale block was not rewritten"
    assert ex.sync_spec_file(path, lemmas) is False, \
        "the fresh block was rewritten again"
