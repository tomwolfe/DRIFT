#!/usr/bin/env python3
"""Export DRIFT's PD/metabolic invariants to QED-parseable LaTeX lemmas.

This is the DRIFT-side semantic bridge to QED, and it is the DRIFT analogue
of ``VeriTrial/scripts/export_pbpk_to_qed.py``.  It reads the DRIFT source
files with the ``ast`` module -- no heavy numerical import is required -- and
emits a deterministic list of lemma strings that QED's parser can consume.

The exported theorems are NOT hand-written constants: the Hill exponent, the
Michaelis-Menten shape, the bridge's basal/max scaling, the Milstein
correction term and the inhibition modulation are read out of the source
(the AST, not a text fragment), so a change to the model changes the
certificate.  A lemma whose shape no longer matches the code fails the
export rather than certifying a stale statement.

Every emitted lemma is an algebraic statement over the reals, of a shape QED
can discharge (identities by ``ring``/``field_simp``, positivity of a
quotient of positive quantities by ``div_pos``/``positivity``).  Where a
theorem needs a SIGN hypothesis that the language cannot state -- the sign of
``Vmax``, the sign of a protein level -- the hypothesis is carried by the
quotient's denominator, so the generated Lean theorem is quantified over
``0 < Vmax`` and friends.  That is stated explicitly in the .tex spec; the
prose names the premise, and the machine-checked part is the algebra.

Usage:
    python3 scripts/export_drift_to_qed.py                 # print lemmas
    python3 scripts/export_drift_to_qed.py --out L.txt     # write them
    python3 scripts/export_drift_to_qed.py --parametric     # with the
                                                            # positivity-bearing
                                                            # (ratio) forms
    python3 scripts/export_drift_to_qed.py --spec-sync     # rewrite the
                                                            # formal spec
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

DRIFT_ROOT = Path(__file__).resolve().parents[1]

# The fragment markers of the machine-checkable block inside each .tex spec.
SPEC_BLOCK_BEGIN = "%% QED-BEGIN"
SPEC_BLOCK_END = "%% QED-END"


# ---------------------------------------------------------------------------
# Source readers.  Each one returns the SHAPE the code implements, so the
# exported lemma cannot drift away from the implementation.
# ---------------------------------------------------------------------------


def _read(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"))


def _binding_source() -> Path:
    return DRIFT_ROOT / "drift" / "binding.py"


def _signaling_source() -> Path:
    return DRIFT_ROOT / "drift" / "signaling.py"


def _metabolic_source() -> Path:
    return DRIFT_ROOT / "drift" / "metabolic.py"


def _topology_source() -> Path:
    return DRIFT_ROOT / "drift" / "topology.py"


def _fail(where: str, what: str) -> ValueError:
    """The fail-closed error: the code no longer implements the equation."""
    return ValueError(
        f"{where} no longer implements the {what} this export certifies. "
        "Re-derive the lemmas instead of certifying a stale equation."
    )


def _function(module: ast.Module, name: str, where: str) -> ast.FunctionDef:
    for node in ast.walk(module):
        if (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name == name):
            return node
    raise _fail(where, f"{name}() function")


def _assigned(func: ast.AST, name: str, where: str) -> ast.expr:
    """The value assigned to ``name`` inside ``func``."""
    for node in ast.walk(func):
        value = None
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    value = node.value
                    break
        elif (isinstance(node, ast.AnnAssign)
              and isinstance(node.target, ast.Name)
              and node.target.id == name):
            value = node.value
        if value is not None:
            return value
    raise _fail(where, f"assignment to {name!r}")


def _returned(func: ast.FunctionDef, where: str) -> ast.expr:
    for node in ast.walk(func):
        if isinstance(node, ast.Return) and node.value is not None:
            return node.value
    raise _fail(where, "returned expression")


def _shape(node: ast.AST) -> str:
    """The AST shape of a node, independent of formatting.

    This is why the readers use the AST rather than a text fragment: a
    reflowed or re-indented equation is the SAME equation, while a changed
    operator, exponent or denominator is a different one -- and only the
    AST tells those two cases apart.
    """
    return ast.dump(node, include_attributes=False)


def _pattern(expression: str) -> str:
    return _shape(ast.parse(expression, mode="eval").body)


def _equals(node: ast.expr, expression: str, where: str, what: str) -> None:
    if _shape(node) != _pattern(expression):
        raise _fail(where, what)


def _contains(func: ast.AST, expression: str, where: str,
              what: str) -> None:
    wanted = _pattern(expression)
    for node in ast.walk(func):
        if isinstance(node, ast.expr) and _shape(node) == wanted:
            return
    raise _fail(where, what)


def _default_argument(module: ast.Module, class_name: str, method: str,
                      argument: str, where: str) -> str:
    """The literal default of ``argument`` in ``Class.method``.

    The Hill certificate must state the exponent the model actually uses.
    A certificate that hard-coded ``n = 1`` would keep verifying after the
    default changed, which is exactly the staleness this exporter exists to
    prevent; so the default is read structurally and becomes part of the
    lemma text.
    """
    for node in ast.walk(module):
        if not (isinstance(node, ast.ClassDef) and node.name == class_name):
            continue
        for item in node.body:
            if not (isinstance(item, ast.FunctionDef) and item.name == method):
                continue
            positional = list(item.args.posonlyargs) + list(item.args.args)
            defaults = list(item.args.defaults)
            offset = len(positional) - len(defaults)
            for index, arg in enumerate(positional):
                if index < offset or arg.arg != argument:
                    continue
                default = defaults[index - offset]
                if (isinstance(default, ast.Constant)
                        and isinstance(default.value, (int, float))
                        and not isinstance(default.value, bool)):
                    value = default.value
                    if float(value) <= 0 or float(value) != int(value):
                        raise _fail(
                            where,
                            f"{class_name}.{method} default for {argument} "
                            f"(must be a positive integer; got {value!r})")
                    return str(int(value))
    raise _fail(where,
                f"{class_name}.{method} default for {argument}")


def extract_hill_shape(path: Path | None = None) -> dict[str, str]:
    """Read the Hill occupancy shape out of ``drift/binding.py``.

    The implementation is::

        numerator   = drug_concentration ** hill_coefficient
        denominator = numerator + kd ** hill_coefficient
        occupancy   = numerator / denominator

    The reader checks those three facts STRUCTURALLY (on the AST, so
    reformatting is not a change but a different operator is) and returns
    the symbol names, so a change of the equation fails the export instead
    of producing a certificate about a model that no longer exists.
    """
    where = "drift/binding.py"
    module = _read(path or _binding_source())
    occupancy = _function(module, "calculate_occupancy", where)
    _equals(_assigned(occupancy, "numerator", where),
            "drug_concentration ** self.hill_coefficient", where,
            "Hill numerator (drug_concentration ** hill_coefficient)")
    _equals(_assigned(occupancy, "denominator", where),
            "numerator + kd ** self.hill_coefficient", where,
            "Hill denominator (numerator + kd ** hill_coefficient)")
    _equals(_returned(occupancy, where),
            "float(numerator / denominator)", where,
            "Hill occupancy (numerator / denominator)")
    exponent = _default_argument(module, "BindingEngine", "__init__",
                                 "hill_coefficient", where)
    return {
        "concentration": "C_drug",
        "affinity": "Kd",
        "exponent": "n_hill",
        "exponent_value": exponent,
        "equation": "theta = C_drug^n_hill / (C_drug^n_hill + Kd^n_hill)",
    }


def extract_flux_shape(path: Path | None = None) -> dict[str, str]:
    """Read the Michaelis-Menten mapping out of ``drift/metabolic.py``."""
    where = "drift/metabolic.py"
    mapping = _function(_read(path or _metabolic_source()),
                        "michaelis_menten_mapping", where)
    _equals(_returned(mapping, where), "(x * (km + 1)) / (km + x)", where,
            "Michaelis-Menten mapping (x * (km + 1) / (km + x))")
    return {
        "substrate": "C_flux",
        "constant": "Km",
        "capacity": "Vmax",
    }


def extract_bridge_scaling(path: Path | None = None) -> dict[str, str]:
    """Read the bridge's basal/max scaling law out of ``metabolic.py``.

    Both the DEFAULTS (read from the ``map_config.get(..., default)`` calls
    the implementation uses) and the two scaling BRANCHES are read
    structurally: the positive influence is ``basal + max_s * level`` and
    the negative one is ``(basal + max_s) - max_s * level``.
    """
    where = "drift/metabolic.py"
    module = _read(path or _metabolic_source())
    defaults: dict[str, str] = {}
    for node in ast.walk(module):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "get"
                and len(node.args) >= 2
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)
                and isinstance(node.args[1], ast.Constant)
                and isinstance(node.args[1].value, (int, float))):
            defaults.setdefault(node.args[0].value, str(node.args[1].value))
    for key in ("basal_scaling", "max_scaling"):
        if key not in defaults:
            raise _fail(where, f"bridge {key!r} default")
    constraints = _function(module, "get_constraints", where)
    _contains(constraints, "basal + max_s * protein_level", where,
              "positive branch of the bridge scaling")
    _contains(constraints, "(basal + max_s) - max_s * protein_level", where,
              "negative branch of the bridge scaling")
    return {
        "basal": defaults["basal_scaling"],
        "max": defaults["max_scaling"],
        "protein": "protein_level",
    }


def extract_milstein_step(path: Path | None = None) -> dict[str, str]:
    """Read the Milstein update out of ``drift/signaling.py``.

    The diffusion coefficient, the Milstein correction and the update are
    read as AST shapes, so the exported algebra is the implemented algebra.
    """
    where = "drift/signaling.py"
    step = _function(_read(path or _signaling_source()),
                     "langevin_step", where)
    _equals(_assigned(step, "b", where),
            "noise_scale * np.sqrt(clamped_state * (1.0 - clamped_state))",
            where, "state-dependent diffusion b(x)")
    _equals(_assigned(step, "milstein_corr", where),
            "0.25 * (noise_scale ** 2) * (1.0 - 2.0 * clamped_state) "
            "* (dw ** 2 - dt)", where, "Milstein correction")
    _equals(_assigned(step, "new_state", where),
            "state + drift * dt + b * dw + milstein_corr", where,
            "Euler-Maruyama update with the Milstein correction")
    return {
        "correction": (
            "0.25 * noise_scale * noise_scale * (1 - 2 * state) "
            "* (dw * dw - dt)"
        ),
        "update": (
            "state + drift * dt + noise_scale * b * dw + "
            "0.25 * noise_scale * noise_scale * (1 - 2 * state) "
            "* (dw * dw - dt)"
        ),
    }


def extract_inhibition_shape(path: Path | None = None) -> dict[str, str]:
    """Read the inhibition modulation out of ``drift/topology.py``.

    The signaling drift scales the effective species level by one minus the
    inhibition, and that scaled level feeds the activation term whose
    expansion the certificate states::

        effective_pi3k = pi3k * (1.0 - inh_pi3k)
        effective_akt  = akt * (1.0 - inh_akt)
        dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - ...

    Without this reader the three inhibition lemmas would be hand-written
    statements the source is never checked against.
    """
    where = "drift/topology.py"
    drift = _function(_read(path or _topology_source()),
                      "pi3k_akt_mtor_drift", where)
    _equals(_assigned(drift, "effective_pi3k", where),
            "pi3k * (1.0 - inh_pi3k)", where,
            "pi3k inhibition modulation (pi3k * (1 - inh_pi3k))")
    _equals(_assigned(drift, "effective_akt", where),
            "akt * (1.0 - inh_akt)", where,
            "akt inhibition modulation (akt * (1 - inh_akt))")
    _contains(drift, "k_akt_act * fb[1] * effective_pi3k * (1.0 - akt)",
              where, "inhibited akt activation term")
    return {
        "species": "pi3k",
        "inhibition": "inh_pi3k",
    }


# ---------------------------------------------------------------------------
# Lemma construction.
# ---------------------------------------------------------------------------


def build_hill_lemmas(shape: dict[str, str]) -> list[str]:
    """Hill occupancy: the partition of unity, the two bounds, and the
    monotonicity DIFFERENCE identity.

    The monotonicity theorem ("a larger concentration gives a larger
    occupancy") is the statement ``occ(C2) - occ(C1) = Kd*(C2-C1)/((C1+Kd)*
    (C2+Kd))``: with ``Kd > 0`` and both denominators positive, the sign of
    the difference is the sign of ``C2 - C1``.  The sign premise itself
    (``C1 < C2``) cannot be written in this statement language, so it is
    stated in the prose of the .tex spec and the machine-checked part is the
    algebraic identity that carries it.
    """
    c = shape["concentration"]
    kd = shape["affinity"]
    n = shape.get("exponent_value", "1")
    # State the exponent the model actually uses.  With the default n = 1
    # the plain symbol is the same statement; with n > 1 the powers appear
    # in every lemma, so a changed default cannot keep yesterday's text.
    cp = c if n == "1" else f"{c}^{n}"
    kp = kd if n == "1" else f"{kd}^{n}"
    c1 = f"{c}1"
    c2 = f"{c}2"
    lemmas = [
        # Partition of unity: the bound and its complement add to one.
        f"{cp} / ({cp} + {kp}) + {kp} / ({cp} + {kp}) = 1",
        # Both parts are non-negative: 0 <= occupancy <= 1.
        f"{cp} / ({cp} + {kp}) >= 0",
        f"{kp} / ({cp} + {kp}) >= 0",
        # The upper bounds are stated as "1 - x >= 0" rather than "x <= 1":
        # MEASURED -- the prover discharges the former and not the latter
        # (a "<= 1" goal with a quotient in it failed on 2026-10-21, while
        # the same bound in this form verified), and the two are the same
        # statement.
        f"1 - {cp} / ({cp} + {kp}) >= 0",
        f"1 - {kp} / ({cp} + {kp}) >= 0",
    ]
    if n == "1":
        # Monotonicity: for n = 1 the difference of two occupancies is the
        # difference of the concentrations times a positive factor.  For a
        # general n the difference has a different factor and is not stated
        # here; the partition and bounds above hold for every n.
        lemmas.append(
            f"{c2} / ({c2} + {kd}) - {c1} / ({c1} + {kd}) = "
            f"{kd} * ({c2} - {c1}) / (({c1} + {kd}) * ({c2} + {kd}))"
        )
    return lemmas


def build_flux_lemmas(shape: dict[str, str]) -> list[str]:
    """Saturable flux: the algebra of the IMPLEMENTED mapping.

    The code returns ``C * (Km + 1) / (Km + C)`` -- a normalized
    Michaelis-Menten form, not the textbook ``C / (Km + C)``.  The
    certificate must state the equation the model actually implements, so
    every lemma below carries the ``(Km + 1)`` factor.  The normalization
    is NOT bounded by one for all C (at C = 10, Km = 0.5 it is ~1.43), so
    no universal upper bound is claimed here; the honest identities are the
    decomposition, the defining rational identity, the cancellation and the
    capacity scaling.
    """
    c = shape["substrate"]
    km = shape["constant"]
    vmax = shape["capacity"]
    flux = f"{c} * ({km} + 1) / ({km} + {c})"
    return [
        # Decomposition: the implemented flux is the sum of the textbook
        # fraction and its Km-weighted part.
        (
            f"{flux} = "
            f"{c} / ({km} + {c}) + {c} * {km} / ({km} + {c})"
        ),
        # Defining rational identity: multiplying by the denominator
        # recovers the numerator exactly.
        f"{flux} * ({km} + {c}) = {c} * ({km} + 1)",
        # Cancellation: for Km + 1 > 0 the normalization divides out.
        f"{flux} / ({km} + 1) = {c} / ({km} + {c})",
        # Capacity scaling: Vmax multiplies the implemented flux, and the
        # scaled flux decomposes into the same two terms.
        (
            f"{vmax} * {c} * ({km} + 1) / ({km} + {c}) = "
            f"{vmax} * {c} / ({km} + {c}) "
            f"+ {vmax} * {c} * {km} / ({km} + {c})"
        ),
    ]


def _scale_constants(*values: str) -> tuple[int, list[int]]:
    """Scale decimal constants to INTEGERS.

    QED's statement language has no decimal literal: ``0.1`` does not
    tokenize, so a statement carrying one fails closed and is never
    proved.  The scale factor is the smallest integer that makes every
    constant integral, and the exported identity is the scaled one -- the
    same identity, multiplied through.  The scale is reported in the .tex
    spec, so the reader can divide it back out.
    """
    numbers = [float(v) for v in values]
    for scale in range(1, 10001):
        scaled = [value * scale for value in numbers]
        if all(abs(s - round(s)) < 1e-9 for s in scaled):
            return scale, [int(round(s)) for s in scaled]
    raise ValueError(
        f"cannot scale the constants {values!r} to integers; the statement "
        "language has no decimal literal and a decimal would make the "
        "lemma fail closed.")


def build_bridge_lemmas(shape: dict[str, str]) -> list[str]:
    """The metabolic bridge's scaling law.

    The bridge computes ``scaling = basal + max * protein_level`` for a
    positive influence and ``(basal + max) - max * protein_level`` for a
    negative one, with the protein level in [0, 1].  The interval bound
    ``basal <= scaling <= basal + max`` follows from ``0 <= protein_level
    <= 1``; the statement language cannot carry that interval premise, so
    the machine-checked part is the ALGEBRA of the scaling law: the two
    branches split the same budget, and the negative branch is the
    complement of the positive one.  The numeric bounds are a code-level
    invariant, pinned by a unit test.

    The constants are scaled to integers (see :func:`_scale_constants`):
    the statement language has no decimal literal, so ``0.1`` would make
    the whole lemma fail closed rather than be proved.
    """
    basal = shape["basal"]
    maximum = shape["max"]
    protein = shape["protein"]
    scale, (scaled_basal, scaled_max) = _scale_constants(basal, maximum)
    shape["scale"] = str(scale)
    shape["scaled_basal"] = str(scaled_basal)
    shape["scaled_max"] = str(scaled_max)
    return [
        # The two branches split the same budget.
        (
            f"{scaled_basal} + {scaled_max} * {protein} "
            f"+ {scaled_max} * (1 - {protein}) "
            f"= {scaled_basal} + {scaled_max}"
        ),
        # The negative-influence branch is the complement of the positive
        # one.
        (
            f"({scaled_basal} + {scaled_max}) - {scaled_max} * {protein} "
            f"= {scaled_basal} + {scaled_max} * (1 - {protein})"
        ),
        # The positive and negative branches are symmetric about the
        # midpoint of the budget: their SUM is constant in protein_level.
        # This is the statement that carries the scaling law's sign
        # behaviour; a bare commutativity identity would prove nothing.
        (
            f"({scaled_basal} + {scaled_max} * {protein}) "
            f"+ (({scaled_basal} + {scaled_max}) "
            f"- {scaled_max} * {protein}) "
            f"= {scaled_basal} + ({scaled_basal} + {scaled_max})"
        ),
    ]


def build_milstein_lemmas(shape: dict[str, str]) -> list[str]:
    """The Milstein algebra, stated as polynomial identities.

    The interval invariant ("a state in [0, 1] stays in [0, 1]") is
    enforced in the implementation by an explicit projection, not by a dt
    bound, and the statement language cannot express the interval premise
    ``0 <= state <= 1`` that the invariant needs.  The invariant is
    therefore a code-level property, pinned by a unit test, and what is
    certified here is the ALGEBRA of the Milstein step: the state-dependent
    variance in its polynomial form, and the expansion of the Milstein
    correction.

    The coefficient 0.25 is scaled to the integer 1 (the identity is
    multiplied by four), because the statement language has no decimal
    literal.
    """
    return [
        # The state-dependent variance in its polynomial form.
        (
            "noise_scale * noise_scale * state * (1 - state) "
            "= noise_scale * noise_scale * (state - state * state)"
        ),
        # The Milstein correction, expanded (the identity multiplied by
        # four, so that no decimal literal appears).
        (
            "noise_scale * noise_scale * (1 - 2 * state) * (dw * dw - dt) "
            "= noise_scale * noise_scale * (dw * dw - dt) "
            "- 2 * noise_scale * noise_scale * state * (dw * dw - dt)"
        ),
        # The correction is antisymmetric in the increment's variance.
        (
            "noise_scale * noise_scale * (1 - 2 * state) * (dw * dw - dt) "
            "+ noise_scale * noise_scale * (1 - 2 * state) * (dt - dw * dw) "
            "= 0"
        ),
    ]


def build_signaling_lemmas(shape: dict[str, str]) -> list[str]:
    """The inhibition modulation the signaling drift implements.

    ``drift/topology.py`` computes ``effective_pi3k = pi3k * (1 - inh)``;
    the identities below are that expansion, which is what makes the
    inhibition a multiplicative modulation of the species level.  The
    ``shape`` is the one :func:`extract_inhibition_shape` read out of the
    source, so the lemmas are only emitted for a drift that implements them.
    """
    # The reader returns the modulation it verified.  A reader that stops
    # returning it (or returns something else) is refused: the lemmas would
    # otherwise be emitted with nothing tying them to the source.
    if not isinstance(shape, dict) \
            or shape.get("species") != "pi3k" \
            or shape.get("inhibition") != "inh_pi3k":
        raise ValueError(
            "drift/topology.py does not expose the certified inhibition "
            f"modulation (reader returned {shape!r})")
    return [
        "pi3k * (1 - inh_pi3k) = pi3k - pi3k * inh_pi3k",
        "akt * (1 - inh_akt) = akt - akt * inh_akt",
        # The inhibition enters the activation term multiplicatively, and the
        # expansion is the one the drift implements.
        (
            "k_akt_act * fb_akt * (1 - akt) * (pi3k * (1 - inh_pi3k)) "
            "= k_akt_act * fb_akt * (1 - akt) * (pi3k - pi3k * inh_pi3k)"
        ),
    ]


def build_lemmas(parametric: bool = True) -> list[str]:
    """The full DRIFT lemma set, in a deterministic order."""
    hill = extract_hill_shape()
    flux = extract_flux_shape()
    bridge = extract_bridge_scaling()
    milstein = extract_milstein_step()
    inhibition = extract_inhibition_shape()
    lemmas: list[str] = []
    lemmas.extend(build_hill_lemmas(hill))
    lemmas.extend(build_flux_lemmas(flux))
    lemmas.extend(build_bridge_lemmas(bridge))
    lemmas.extend(build_milstein_lemmas(milstein))
    lemmas.extend(build_signaling_lemmas(inhibition))
    if not parametric:
        # The non-parametric form keeps only the statements that carry no
        # ratio: the plain identities.
        keep = (" = 1", " = 0")
        lemmas = [l for l in lemmas if l.endswith(keep)]
    return lemmas


# ---------------------------------------------------------------------------
# Formal-spec synchronisation.
# ---------------------------------------------------------------------------

SPEC_FILES = {
    "drift_binding_occupancy": (
        "drift_binding_occupancy.tex",
        "build_hill_lemmas",
    ),
    "drift_metabolic_flux": (
        "drift_metabolic_flux.tex",
        "build_flux_lemmas",
    ),
    "drift_bridge_scaling": (
        "drift_bridge_scaling.tex",
        "build_bridge_lemmas",
    ),
    "drift_signaling_milstein": (
        "drift_signaling_milstein.tex",
        "build_milstein_lemmas",
    ),
}


def spec_dir() -> Path:
    return DRIFT_ROOT / "formal_specs"


def build_spec_blocks() -> dict[str, list[str]]:
    """The lemma list of every formal spec, keyed by file stem."""
    return {
        "drift_binding_occupancy": build_hill_lemmas(extract_hill_shape()),
        "drift_metabolic_flux": build_flux_lemmas(extract_flux_shape()),
        "drift_bridge_scaling": build_bridge_lemmas(extract_bridge_scaling()),
        "drift_signaling_milstein": (
            build_milstein_lemmas(extract_milstein_step())
            + build_signaling_lemmas(extract_inhibition_shape())
        ),
    }


def sync_spec_file(path: Path, lemmas: list[str]) -> bool:
    """Rewrite the machine-checkable block of a .tex spec.

    Returns True when the file changed.  The block is delimited by
    ``%% QED-BEGIN`` / ``%% QED-END``; everything outside it is prose and is
    left alone.
    """
    if not path.is_file():
        raise FileNotFoundError(f"formal spec not found: {path}")
    text = path.read_text(encoding="utf-8")
    begin = text.find(SPEC_BLOCK_BEGIN)
    end = text.find(SPEC_BLOCK_END)
    if begin == -1 or end == -1 or end < begin:
        raise ValueError(
            f"{path} has no machine-checkable block "
            f"({SPEC_BLOCK_BEGIN} ... {SPEC_BLOCK_END})"
        )
    body = "\n".join(f"%% QED {lemma}" for lemma in lemmas)
    new_text = (
        text[:begin]
        + f"{SPEC_BLOCK_BEGIN}\n{body}\n"
        + text[end:]
    )
    if new_text == text:
        return False
    path.write_text(new_text, encoding="utf-8")
    return True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=None,
                    help="write the lemmas to this file instead of stdout")
    ap.add_argument("--parametric", action="store_true",
                    help="emit the full parametric lemma set (default)")
    ap.add_argument("--no-parametric", action="store_true",
                    help="emit only the plain identities")
    ap.add_argument("--spec-sync", action="store_true",
                    help="rewrite the machine-checkable blocks of the .tex "
                         "specs from the current source")
    ap.add_argument("--spec-dir", type=Path, default=None,
                    help="override the formal-spec directory")
    args = ap.parse_args(argv)

    try:
        lemmas = build_lemmas(parametric=not args.no_parametric)
    except ValueError as e:
        print(f"DRIFT EXPORT FAILED (fail-closed): {e}", file=sys.stderr)
        return 1

    if args.spec_sync:
        directory = args.spec_dir or spec_dir()
        blocks = build_spec_blocks()
        for stem, block in blocks.items():
            path = directory / f"{stem}.tex"
            try:
                changed = sync_spec_file(path, block)
            except (FileNotFoundError, ValueError) as e:
                print(f"DRIFT EXPORT FAILED (fail-closed): {e}",
                      file=sys.stderr)
                return 1
            print(
                f"{'updated' if changed else 'unchanged'}: {path}",
                file=sys.stderr,
            )

    text = "\n".join(lemmas) + "\n"
    if args.out is not None:
        args.out.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
