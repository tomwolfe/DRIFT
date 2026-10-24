#!/usr/bin/env python3
"""End-to-end in-silico trial: VeriTrial PK -> DRIFT PD/metabolism.

This is the seam that turns the two engines into ONE trial:

    VeriTrial PBPK  --(C_p(t))-->  DRIFT Hill binding
                     --(inhibition(t))-->  DRIFT stochastic signaling
                     --(state)-->  DRIFT metabolic bridge
                     --(constraints)-->  DRIFT FBA  -->  growth / cell death

VeriTrial is the SAD/MAD engine (it owns the concentration-time profile and
the regulatory NCA/PD summaries); DRIFT is the PD/metabolic engine (it owns
the target occupancy, the signaling network, the metabolic constraints and
the growth readout).  The seam between them is the PLASMA CONCENTRATION:
the engine is driven by VeriTrial's C_p(t) profile, one concentration per
step, and never by a hand-fed scalar.

What the script ASSERTS, and how -- every one of them is MEASURED:

  A. the inhibition is monotone in the concentration, measured on the real
     C_p(t) profile of the trial dose: the sorted profile has a sorted
     inhibition, and the SIGN of the occupancy difference is the sign of
     the concentration difference -- the sign reasoning of the algebra the
     formal layer certifies (the Hill difference identity);
  B. ONE trial: the engine's own recorded concentration series IS the
     plasma profile, and its own recorded inhibition series IS the Hill
     occupancy of that profile;
  C. the growth readout is monotone in the DOSE (in mg), measured end to
     end through VeriTrial's PBPK and DRIFT's engine: more drug means less
     growth;
  D. the cell_death flag is reachable, and is reached by driving a genuinely
     INFEASIBLE metabolic problem -- the LP solver's own verdict, observed
     by the engine, never a flag set by hand.  A feasible control with the
     same wiring must NOT die, so the flag tracks the infeasibility;
  E. the formal layer still verifies through QED in this same run: the Hill
     bounds and monotonicity algebra and the IMPLEMENTED normalized flux
     identities, taken from DRIFT's own exporter, vetted with the DRIFT
     gate's own rules and PROVED by QED (no sorry, nothing skipped); the
     flux family is not called "bounds" because C*(Km+1)/(Km+C) is not
     bounded by one;
  F. the provenance file chains the FOUR repository SHAs, the Lean artefact
     digests and a benchmark summary.

The script FAILS CLOSED: any failed assertion exits non-zero, and the
artefacts are written only once every assertion has passed.

Usage:
    python3 scripts/run_veritrial_drift_integration.py
    python3 scripts/run_veritrial_drift_integration.py --out output/integration.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

import numpy as np

DRIFT_ROOT = Path(__file__).resolve().parents[1]
APPS_ROOT = DRIFT_ROOT.parent
VERITRIAL_ROOT = APPS_ROOT / "VeriTrial"
QED_ROOT = APPS_ROOT / "QED"
TETHER_ROOT = APPS_ROOT / "tether"

# The trial: one 400 mg oral moxifloxacin dose, observed for 24 h.
TRIAL_DOSE_MG = 400.0
TRIAL_HOURS = 24.0
PROFILE_POINTS = 49

# The engine's own time step (h) and the number of steps that cover the
# trial: 24 h / 0.05 h = 480 steps.
ENGINE_DT = 0.05
TRIAL_STEPS = int(round(TRIAL_HOURS / ENGINE_DT))

# The dose ladder the growth monotonicity is measured on, in mg.
DOSE_LADDER_MG = (0.0, 100.0, 200.0, 400.0, 800.0)

# One seed for every run of the ladder: the doses are compared PAIRWISE,
# so the difference between two doses is the drug, not the noise.
SEED = 20261021

# The biomass DEMAND the cell-death run drives.  The e_coli_core network
# can deliver at most ~0.87 mmol/gDW/h of biomass, so a demand of 5.0 is
# infeasible; the control's demand of 0.1 is comfortably feasible.
INFEASIBLE_DEMAND = 5.0
FEASIBLE_DEMAND = 0.1
BIOMASS_REACTION = "Biomass_Ecoli_core"

# The Lean artefacts the chain digests.  ``output.lean`` is deliberately
# NOT among them: QED rewrites it on every proof, so its digest is
# volatile and would make the chain unstable.
LEAN_ARTEFACTS = (
    "Compartmental.lean",
    "QED.lean",
    "VeriTrialExport.lean",
    "lakefile.lean",
)

THE_FOUR_REPOS = ("DRIFT", "VeriTrial", "QED", "tether")


class IntegrationFailure(RuntimeError):
    """One failed assertion.  The run fails closed on it."""


def require(condition: bool, message: str) -> None:
    """Assert one measured property of the trial."""
    if not condition:
        raise IntegrationFailure(message)


def _add_to_path(path: Path) -> None:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


# ---------------------------------------------------------------------------
# VeriTrial: the concentration-time profile.
# ---------------------------------------------------------------------------


def central_concentration_profile(
    dose_mg: float = TRIAL_DOSE_MG,
    hours: float = TRIAL_HOURS,
    points: int = PROFILE_POINTS,
    weight_kg: float = 70.0,
    age: int = 30,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """VeriTrial's central-compartment concentration for one dose.

    Returns ``(t_hours, C_p, params)``.  The absorbed amount is the
    administered dose scaled by the drug's bioavailability, which is what
    VeriTrial's own patient-level path does.
    """
    _add_to_path(VERITRIAL_ROOT / "src")
    from insilico_trial.pbpk.model import build_pbpk_params, solve_pbpk_single
    from insilico_trial.schemas import load_drug_config

    drug = load_drug_config(VERITRIAL_ROOT / "configs" / "drug_moxifloxacin.yaml")
    params = build_pbpk_params(
        weight_kg=weight_kg,
        age=age,
        drug=drug,
        genotype_scale=1.0,
    )
    t_eval = np.linspace(0.0, hours, points)
    absorbed = dose_mg * float(drug.bioavailability)
    C_p = np.asarray(solve_pbpk_single(t_eval, absorbed, params), dtype=np.float64)
    return t_eval, C_p, {"params": params, "drug": drug}


def read_drug_parameters() -> dict:
    """The drug's own PD parameters, read from VeriTrial's configuration."""
    import yaml  # noqa: PLC0415

    drug_cfg = yaml.safe_load(
        (VERITRIAL_ROOT / "configs" / "drug_moxifloxacin.yaml").read_text())
    drug = drug_cfg["drug"]
    return {
        "ec50": float(drug["ec50"]),
        "hill_coeff": float(drug["hill_coeff"]),
        "bioavailability": float(drug["bioavailability"]),
    }


# ---------------------------------------------------------------------------
# The Hill occupancy, written out here so the measurement does not depend
# on the engine's own arithmetic.
# ---------------------------------------------------------------------------


def hill_occupancy(concentration: float, kd: float, hill: float) -> float:
    """``theta = C^n / (C^n + Kd^n)``, the occupancy DRIFT's engine uses."""
    if concentration <= 0.0:
        return 0.0
    numerator = concentration ** hill
    denominator = numerator + kd ** hill
    return float(numerator / denominator)


def hill_difference(c1: float, c2: float, kd: float, n: float) -> float:
    """``occ(c2) - occ(c1)`` in the closed form the formal layer certifies.

    The certified lemma is the ``n = 1`` instance::

        occ(c2) - occ(c1) = Kd * (c2 - c1) / ((c1 + Kd) * (c2 + Kd))

    and the model's Hill coefficient ``n`` gives the same identity with
    every power raised to ``n``::

        occ(c2) - occ(c1) = Kd^n (c2^n - c1^n) / ((c1^n + Kd^n) (c2^n + Kd^n))

    With ``Kd > 0`` and both denominators positive, the sign of the
    difference is the sign of ``c2 - c1``.
    """
    if n == 1.0:
        return kd * (c2 - c1) / ((c1 + kd) * (c2 + kd))
    return (kd ** n * (c2 ** n - c1 ** n)
            / ((c1 ** n + kd ** n) * (c2 ** n + kd ** n)))


# ---------------------------------------------------------------------------
# DRIFT: the engine.
# ---------------------------------------------------------------------------


def build_engine(km: float, hill: float, bridge: Any = None,
                 dt: float = ENGINE_DT):
    """DRIFT's PD/metabolic engine, wired to the Hill binding."""
    _add_to_path(DRIFT_ROOT)
    from drift.binding import BindingEngine
    from drift.engine import SimulationEngine
    from drift.metabolic import DFBASolver, MetabolicBridge
    from drift.signaling import StochasticIntegrator
    from drift.topology import get_default_topology

    topology = get_default_topology()
    integrator = StochasticIntegrator(dt=dt, noise_scale=0.01, topology=topology)
    solver = DFBASolver(model_name="textbook")
    binding = BindingEngine(targets=km, hill_coefficient=hill)
    return SimulationEngine(
        integrator, solver,
        bridge if bridge is not None else MetabolicBridge(),
        binding,
    )


def run_trial(km: float, hill: float, dose_mg: float,
              steps: int = TRIAL_STEPS, bridge: Any = None,
              seed: int = SEED) -> dict:
    """Run ONE trial: VeriTrial's PK for the dose, then DRIFT's PD and
    metabolism over the whole profile."""
    t_hours, C_p, _ = central_concentration_profile(
        dose_mg=dose_mg, hours=TRIAL_HOURS)
    engine = build_engine(km, hill, bridge=bridge)
    result = engine.run(steps=steps, drug_profile=C_p,
                         drug_profile_times=t_hours, seed=seed)
    growth = np.asarray(result["growth"], dtype=np.float64)
    times = np.asarray(result["time"], dtype=np.float64)
    return {
        "dose_mg": float(dose_mg),
        "cmax_mg_per_l": float(C_p.max()),
        "conc_auc_mg_h_per_l": float(np.trapezoid(C_p, t_hours)),
        "inhibition_max": float(np.max(result["inhibition_series"])),
        "times": times,
        "profile_times": t_hours,
        "profile": C_p,
        "concentration": np.asarray(result["concentration"], dtype=np.float64),
        "inhibition_series": np.asarray(
            result["inhibition_series"], dtype=np.float64),
        "mean_growth": float(growth.mean()) if growth.size else 0.0,
        "growth_auc": (float(np.trapezoid(growth, times))
                       if growth.size else 0.0),
        "max_growth": float(growth.max()) if growth.size else 0.0,
        "cell_death": bool(result["cell_death"]),
        "death_step": result["death_step"],
        "death_cause": result["death_cause"],
        "statuses": list(result["status"]),
    }


def requirement_bridge(demand: float):
    """A bridge that states a DEMAND on the biomass reaction.

    A capacity mapping caps what a reaction may do; a requirement
    mapping states what a reaction MUST carry.  The demand is the value
    the biomass flux has to reach, so a demand above what the network
    can produce makes the LP INFEASIBLE -- which is how a metabolic
    problem reports that the cell cannot sustain itself.  The glucose
    capacity mapping around it is the default one, so the wiring of the
    demand is the wiring of an ordinary run.
    """
    _add_to_path(DRIFT_ROOT)
    from drift.metabolic import MetabolicBridge

    return MetabolicBridge(
        mappings=[
            {
                "protein_name": "mTOR",
                "reaction_id": "EX_glc__D_e",
                "influence": "positive",
                "base_vmax": 10.0,
                "basal_scaling": 0.1,
                "max_scaling": 0.9,
            },
            {
                "protein_name": "mTOR",
                "reaction_id": BIOMASS_REACTION,
                "constraint_kind": "requirement",
                "influence": "positive",
                "base_vmax": demand,
                "basal_scaling": 1.0,
                "max_scaling": 0.0,
            },
        ],
        species_names=["PI3K", "AKT", "mTOR"],
    )


def drive_metabolic_demand(km: float, hill: float, demand: float) -> dict:
    """Drive a metabolic problem with a biomass demand, and report what
    the LP solver and the engine each made of it.

    Nothing is set by hand: the bridge states the demand, the LP solver
    returns its own verdict, and the engine observes that verdict.
    """
    bridge = requirement_bridge(demand)
    engine = build_engine(km, hill, bridge=bridge)

    # The same problem, solved directly, so the infeasibility is the
    # SOLVER's verdict and not the engine's opinion of it.
    state = engine.integrator.topology.get_initial_state()
    constraints, scalings = bridge.get_constraints_with_scalings(state)
    direct = engine.solver.optimize(constraints, scalings=scalings)

    result = engine.run(steps=4, drug_concentration=1.0, seed=3)
    growth = np.asarray(result["growth"], dtype=np.float64)
    return {
        "demand": float(demand),
        "constraints": {k: float(v) for k, v in constraints.items()},
        "lp_status": direct.get("status"),
        "lp_objective": float(direct.get("objective_value", 0.0)),
        "lp_bottlenecks": dict(direct.get("bottlenecks", {}) or {}),
        "lp_diagnostic": direct.get("diagnostic"),
        "cell_death": bool(result["cell_death"]),
        "death_step": result["death_step"],
        "death_cause": result["death_cause"],
        "growth": [float(g) for g in growth],
        "statuses": list(result["status"]),
    }


# ---------------------------------------------------------------------------
# The formal layer.
# ---------------------------------------------------------------------------


def verify_formal_layer() -> dict:
    """Verify the mission's formal layer through QED, in this same run.

    The three families the mission names -- the Hill bounds and
    monotonicity algebra, and the IMPLEMENTED normalized flux identities --
    are taken from DRIFT's own exporter, vetted with the DRIFT gate's own
    rules (the set must be exactly what the exporter emits, no theater, no
    sorry) and then PROVED by QED.  A lemma that is skipped instead of
    proved fails the gate, exactly as it does in the standalone formal gate.

    The flux family is deliberately NOT called "bounds": the implemented
    mapping is C*(Km+1)/(Km+C), which exceeds one for large C, so claiming
    0 <= flux <= 1 would certify a false statement.  The family certifies
    the algebraic identities the code actually satisfies.
    """
    _add_to_path(DRIFT_ROOT / "scripts")
    import export_drift_to_qed as exporter  # noqa: PLC0415
    import verify_drift_lemmas as gate  # noqa: PLC0415

    families = {
        "hill_bounds_and_monotonicity": exporter.build_hill_lemmas(
            exporter.extract_hill_shape()),
        "flux_normalized_mapping": exporter.build_flux_lemmas(
            exporter.extract_flux_shape()),
    }
    lemmas = [lemma for block in families.values() for lemma in block]

    # The named families must really be in there: Hill bounds that are
    # inequalities and the monotonicity DIFFERENCE identity; flux
    # identities that carry the implemented (Km + 1) normalization.
    def _is_monotonicity(lemma: str) -> bool:
        """The difference identity: two occupancies on the left, a single
        quotient with a difference in its numerator on the right."""
        left, _, right = lemma.partition("=")
        return (lemma.count("=") == 1
                and " - " in left
                and " - " in right
                and right.strip().endswith(")"))

    hill_bounds = [l for l in families["hill_bounds_and_monotonicity"]
                   if l.rstrip().endswith(">= 0")]
    flux_identities = [l for l in families["flux_normalized_mapping"]
                       if l.count("=") == 1 and "(Km + 1)" in l]
    monotonicity = [l for l in families["hill_bounds_and_monotonicity"]
                    if _is_monotonicity(l)]
    require(len(hill_bounds) >= 2,
            "the Hill family carries no occupancy bounds")
    require(len(flux_identities) >= 2,
            "the flux family does not certify the implemented "
            "(Km + 1) mapping")
    require(len(monotonicity) == 1,
            "the Hill family carries no monotonicity difference identity")

    # The gate's own teeth: the set must be what DRIFT emits right now,
    # no lemma may be theater, and no lemma may carry a sorry.
    live = gate.live_lemmas()
    missing = sorted(set(lemmas) - set(live))
    require(not missing,
            f"the mission's formal families are not what the exporter "
            f"emits: {missing}")
    theater = [l for l in lemmas if gate._is_theater(l)]
    require(not theater,
            f"the lemma set contains statements that prove nothing about "
            f"the model: {theater}")
    sorry = [l for l in lemmas if gate._is_sorry_placeholder(l)]
    require(not sorry,
            f"the lemma set carries a sorry placeholder: {sorry}")

    verifier = QED_ROOT / "verify_pbpk_lemmas.py"
    require(verifier.is_file(),
            f"QED not found at {QED_ROOT}; the formal layer cannot be "
            "verified, so the gate fails closed")

    with tempfile.NamedTemporaryFile("w", suffix=".lemmas",
                                     delete=False) as handle:
        handle.write("\n".join(lemmas) + "\n")
        lemma_file = handle.name
    try:
        proc = subprocess.run(
            [sys.executable, str(verifier), lemma_file],
            cwd=str(QED_ROOT), capture_output=True, text=True, timeout=3600,
        )
    finally:
        Path(lemma_file).unlink(missing_ok=True)

    require(proc.returncode == 0,
            f"QED did not verify the mission's formal layer:\n"
            f"{proc.stdout}\n{proc.stderr}")
    for lemma in lemmas:
        require(f"verified: {lemma}" in proc.stdout,
                f"QED did not verify: {lemma}")
    skipped = [line for line in proc.stdout.splitlines()
               if line.startswith("skipped")]
    require(not skipped,
            f"QED skipped a lemma instead of proving it: {skipped}")
    require("no sorry" in proc.stdout,
            "QED's own summary does not report a sorry-free run")

    return {
        "qed": str(QED_ROOT),
        "families": {name: list(block) for name, block in families.items()},
        "lemmas": lemmas,
        "verified": len(lemmas),
        "skipped": 0,
        "sorry": 0,
        "output": [line for line in proc.stdout.splitlines()
                   if line.startswith("verified:")],
    }


# ---------------------------------------------------------------------------
# The provenance chain.
# ---------------------------------------------------------------------------


def _git_sha(repo: Path) -> Optional[str]:
    proc = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return None
    return proc.stdout.strip() or None


def _git_dirty(repo: Path) -> Optional[bool]:
    proc = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return None
    return bool(proc.stdout.strip())


def _sha256(path: Path) -> str:
    if not path.is_file():
        return "missing"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_tether_ledger() -> dict:
    """The tether ledger: the cross-repo record of the four heads."""
    ledger = TETHER_ROOT / "SYSTEM_STATE.json"
    if not ledger.is_file():
        return {}
    try:
        data = json.loads(ledger.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def drift_source_files() -> list[Path]:
    """The DRIFT sources that produced this run."""
    files = sorted((DRIFT_ROOT / "drift").rglob("*.py"))
    files += sorted((DRIFT_ROOT / "scripts").glob("*.py"))
    return [p for p in files if "__pycache__" not in p.parts]


def drift_source_digest() -> str:
    """A content digest of the DRIFT sources that ran.

    It binds a provenance record to the code that produced it, so a
    record carried over from another environment is only accepted when
    it is about the same sources.
    """
    digest = hashlib.sha256()
    for path in drift_source_files():
        digest.update(path.relative_to(DRIFT_ROOT).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _ledger_head(ledger: dict, name: str) -> Optional[str]:
    for row in ledger.get("repos", []) or []:
        if str(row.get("repo")) == name:
            head = row.get("head")
            return str(head) if head else None
    return None


def resolve_repo_head(name: str, path: Path, ledger: dict,
                      carried: dict) -> dict:
    """The head of one repository, and where the answer came from.

    The chain names four repositories.  A repository that is not a git
    checkout in this environment -- the clean room carries the project
    itself as an archive, without its ``.git`` -- is resolved from the
    tether ledger when the ledger is reachable, and otherwise from the
    record the run that produced these sources left behind.  That
    carried record is accepted for DRIFT only when it carries the digest
    of the sources that are running now.
    """
    sha = _git_sha(path)
    if sha:
        return {"sha": sha, "sha_source": "git"}
    head = _ledger_head(ledger, name)
    if head:
        return {"sha": head, "sha_source": "tether-ledger"}
    prior = (carried.get("repos") or {}).get(name) or {}
    prior_sha = prior.get("sha")
    if prior_sha and prior_sha != "unknown":
        binding = str(carried.get("drift_source_digest") or "")
        if name == "DRIFT" and binding != drift_source_digest():
            # The carried record is about other sources: it does not
            # describe this run.
            return {"sha": "unknown", "sha_source": "unresolved"}
        return {
            "sha": prior_sha,
            "sha_source": "carried-forward",
            "carried_from": "output/regulatory_provenance.json",
            "binding": binding,
        }
    return {"sha": "unknown", "sha_source": "unresolved"}


def build_provenance(payload: dict, formal: dict, carried: dict) -> dict:
    """The provenance chain: four repository SHAs, the Lean artefact
    digests and a benchmark summary."""
    ledger = load_tether_ledger()
    repos = {
        "DRIFT": DRIFT_ROOT,
        "VeriTrial": VERITRIAL_ROOT,
        "QED": QED_ROOT,
        "tether": TETHER_ROOT,
    }
    repo_rows: dict[str, dict] = {}
    for name, path in repos.items():
        row = resolve_repo_head(name, path, ledger, carried)
        row["dirty"] = _git_dirty(path)
        repo_rows[name] = row

    lean_digests = {
        f"QED/{name}": _sha256(QED_ROOT / name)
        for name in LEAN_ARTEFACTS
    }
    benchmark = payload.get("benchmark", {})

    chain = {
        "repos": {name: row["sha"] for name, row in repo_rows.items()},
        "lean_digests": lean_digests,
        "benchmark": benchmark,
    }
    merkle_root = hashlib.sha256(
        json.dumps(chain, sort_keys=True).encode("utf-8")).hexdigest()

    # Containment with the tether ledger: the same four heads and the
    # same Lean digests appear in both.  The ledger hashes a DIFFERENT
    # projection (its own audit rows), so the two roots are NOT expected
    # to be equal -- asserting equality would be the bug.
    containment: dict[str, Any] = {"ledger_available": bool(ledger)}
    if ledger:
        heads_match = {}
        for name in THE_FOUR_REPOS:
            ledger_head = _ledger_head(ledger, name)
            heads_match[name] = (
                ledger_head is not None
                and ledger_head == repo_rows[name]["sha"]
            )
        digests_match = {}
        for row in ledger.get("repos", []) or []:
            if str(row.get("repo")) != "QED":
                continue
            for name, digest in (row.get("lean_files") or {}).items():
                if name == "output.lean":
                    # Volatile: QED rewrites it on every proof.
                    continue
                key = f"QED/{name}"
                if key in lean_digests:
                    digests_match[key] = lean_digests[key] == digest
        containment.update({
            "heads_match": heads_match,
            "lean_digests_match": digests_match,
            "ledger_merkle_root": ledger.get("merkle_root"),
        })
        require(all(heads_match.values()),
                f"the provenance heads do not match the tether ledger: "
                f"{heads_match}")
        require(all(digests_match.values()),
                f"the Lean digests do not match the tether ledger: "
                f"{digests_match}")
        require(bool(ledger.get("merkle_root")),
                "the tether ledger carries no merkle root")

    return {
        "schema": "drift-veritrial-integration-provenance/2",
        "merkle_root": merkle_root,
        "tether_ledger_merkle_root": ledger.get("merkle_root"),
        "tether_ledger": containment,
        "repos": repo_rows,
        "lean_digests": lean_digests,
        "benchmark": benchmark,
        "drift_source_digest": drift_source_digest(),
        "formal": {
            "qed": formal.get("qed"),
            "verified": formal.get("verified"),
            "skipped": formal.get("skipped"),
            "sorry": formal.get("sorry"),
            "families": formal.get("families"),
        },
        "integration": {
            key: payload[key]
            for key in ("dose_mg", "km", "hill", "cmax", "pmax",
                        "inhibition_max", "dose_ladder", "cell_death",
                        "formal_gate")
            if key in payload
        },
    }


def load_carried_provenance(path: Path) -> dict:
    """The provenance record a previous run left, if there is one."""
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


# ---------------------------------------------------------------------------
# The assertions.
# ---------------------------------------------------------------------------


def assert_inhibition_is_monotone(C_p: np.ndarray, km: float,
                                   hill: float) -> dict:
    """A. The inhibition is monotone in the concentration, MEASURED on the
    real C_p(t) profile of the trial dose."""
    require(C_p.size >= 2,
            "the concentration profile carries no samples")
    inhibition = np.asarray(
        [hill_occupancy(float(c), km, hill) for c in C_p], dtype=np.float64)

    order = np.argsort(C_p, kind="stable")
    sorted_inhibition = inhibition[order]
    require(bool(np.all(np.diff(sorted_inhibition) >= -1e-12)),
            "the Hill inhibition is not monotone in the concentration")
    require(0.0 <= float(inhibition.min()) and float(inhibition.max()) <= 1.0,
            "the inhibition left [0, 1]")

    # The sign link to the algebra the formal layer certifies: the
    # occupancy difference FOLLOWS the certified identity, and its sign
    # is the sign of the concentration difference -- the reasoning the
    # certified Hill difference identity carries.  EVERY pair of
    # distinct concentrations on the profile is checked, in both
    # directions.
    pairs = 0
    for i in range(C_p.size):
        for j in range(i + 1, C_p.size):
            c1, c2 = float(C_p[i]), float(C_p[j])
            if c2 == c1:
                continue
            measured = (hill_occupancy(c2, km, hill)
                        - hill_occupancy(c1, km, hill))
            certified = hill_difference(c1, c2, km, hill)
            tolerance = 1e-9 * max(abs(measured), abs(certified), 1e-12)
            require(abs(measured - certified) <= tolerance,
                    f"the occupancy difference at C={c1}, C={c2} does not "
                    f"follow the certified identity ({measured} != "
                    f"{certified})")
            certified_n1 = hill_difference(c1, c2, km, 1.0)
            require(bool(np.sign(measured) == np.sign(c2 - c1)
                         == np.sign(certified_n1)),
                    f"the sign of the occupancy difference at C={c1}, C={c2} "
                    "is not the sign of the concentration difference")
            pairs += 1
    require(pairs > 0,
            "the concentration profile carries no two distinct "
            "concentrations, so monotonicity is vacuous")
    return {
        "pairs_checked": pairs,
        "inhibition_min": float(inhibition.min()),
        "inhibition_max": float(inhibition.max()),
        "monotone": True,
    }


def assert_one_trial(trial: dict, t_hours: np.ndarray, C_p: np.ndarray,
                      km: float, hill: float) -> dict:
    """B. ONE trial: the engine's own record IS the plasma profile and the
    Hill occupancy of that profile."""
    expected_c = np.interp(trial["times"], t_hours, C_p)
    require(bool(np.allclose(trial["concentration"], expected_c,
                              rtol=0.0, atol=1e-12)),
            "the engine did not see the plasma concentration profile")
    expected_inhibition = np.asarray(
        [hill_occupancy(float(c), km, hill) for c in expected_c],
        dtype=np.float64)
    require(bool(np.allclose(trial["inhibition_series"], expected_inhibition,
                              rtol=0.0, atol=1e-12)),
            "the engine's inhibition is not the Hill occupancy of the "
            "plasma profile")
    return {
        "steps": int(trial["times"].size),
        "concentration_matches_profile": True,
        "inhibition_matches_hill": True,
    }


def assert_growth_is_monotone_in_dose(km: float, hill: float) -> list[dict]:
    """C. More drug means less growth, measured end to end through the
    engine (in mg, not in a hand-fed concentration)."""
    ladder = [run_trial(km, hill, dose) for dose in DOSE_LADDER_MG]
    for row in ladder:
        require(not row["cell_death"],
                f"the {row['dose_mg']} mg run died: {row['death_cause']}")
    for lower, higher in zip(ladder, ladder[1:]):
        require(higher["mean_growth"] < lower["mean_growth"],
                f"the growth readout is not monotone in the dose: "
                f"{lower['dose_mg']} mg -> {higher['dose_mg']} mg gives "
                f"{lower['mean_growth']} -> {higher['mean_growth']}")
        require(higher["growth_auc"] <= lower["growth_auc"],
                f"the growth AUC is not monotone in the dose: "
                f"{lower['dose_mg']} mg -> {higher['dose_mg']} mg gives "
                f"{lower['growth_auc']} -> {higher['growth_auc']}")
    return ladder


def assert_cell_death_is_reachable(km: float, hill: float) -> dict:
    """D. The cell-death flag is reached by an INFEASIBLE metabolic
    problem, and the flag tracks the infeasibility."""
    death = drive_metabolic_demand(km, hill, INFEASIBLE_DEMAND)
    require(death["lp_status"] == "infeasible",
            f"the metabolic problem was not infeasible: {death['lp_status']}")
    require(BIOMASS_REACTION in death["lp_bottlenecks"],
            f"the infeasibility is not the one that was driven: "
            f"{death['lp_bottlenecks']}")
    require(death["cell_death"],
            "the engine did not set the cell-death flag for an infeasible "
            "metabolic problem")
    require(death["death_step"] == 0,
            f"the cell death was not detected at the first step: "
            f"{death['death_step']}")
    require(all(growth == 0.0 for growth in death["growth"]),
            f"the growth was not zeroed after the death: {death['growth']}")
    require(BIOMASS_REACTION in (death["death_cause"] or ""),
            f"the death cause is not the solver's own diagnostic: "
            f"{death['death_cause']}")

    control = drive_metabolic_demand(km, hill, FEASIBLE_DEMAND)
    require(control["lp_status"] == "optimal",
            f"the control metabolic problem was not feasible: "
            f"{control['lp_status']}")
    require(not control["cell_death"],
            "the cell-death flag does not track the infeasibility: the "
            "feasible control died as well")
    require(control["growth"] and all(growth > 0.0 for growth in control["growth"]),
            f"the feasible control grew nothing: {control['growth']}")
    return {"infeasible": death, "feasible_control": control}


# ---------------------------------------------------------------------------
# The run.
# ---------------------------------------------------------------------------


def run_integration(args: argparse.Namespace) -> tuple[dict, dict]:
    """Every assertion, in order.  Raises on the first failure."""
    drug = read_drug_parameters()
    km = drug["ec50"]
    hill = drug["hill_coeff"]

    # A. The inhibition is monotone in the concentration.
    t_hours, C_p, _ = central_concentration_profile(dose_mg=args.dose)
    monotonicity = assert_inhibition_is_monotone(C_p, km, hill)

    # B. ONE trial: the engine runs on the profile.
    trial = run_trial(km, hill, args.dose)
    one_trial = assert_one_trial(trial, t_hours, C_p, km, hill)

    # C. More drug means less growth.
    ladder = assert_growth_is_monotone_in_dose(km, hill)

    # D. The cell-death flag is reachable through a real infeasibility.
    death = assert_cell_death_is_reachable(km, hill)

    # E. The formal layer, in this same run.
    formal = verify_formal_layer()

    ladder_rows = [
        {
            "dose_mg": row["dose_mg"],
            "cmax_mg_per_l": row["cmax_mg_per_l"],
            "conc_auc_mg_h_per_l": row["conc_auc_mg_h_per_l"],
            "inhibition_max": row["inhibition_max"],
            "mean_growth": row["mean_growth"],
            "growth_auc": row["growth_auc"],
            "max_growth": row["max_growth"],
            "cell_death": row["cell_death"],
        }
        for row in ladder
    ]

    payload = {
        "dose_mg": float(args.dose),
        "km": km,
        "hill": hill,
        "bioavailability": drug["bioavailability"],
        "cmax": float(C_p.max()),
        "pmax": float(C_p.max()),
        "conc_auc": float(np.trapezoid(C_p, t_hours)),
        "inhibition_max": float(np.max(trial["inhibition_series"])),
        "hours": TRIAL_HOURS,
        "steps": TRIAL_STEPS,
        "seed": SEED,
        "monotonicity": monotonicity,
        "one_trial": one_trial,
        "dose_ladder": ladder_rows,
        "cell_death": death,
        "formal_gate": {
            "qed": formal["qed"],
            "verified": formal["verified"],
            "skipped": formal["skipped"],
            "sorry": formal["sorry"],
            "lemmas": formal["lemmas"],
        },
        "benchmark": {
            "growth_by_dose": {str(row["dose_mg"]): row["mean_growth"]
                               for row in ladder_rows},
            "growth_auc_by_dose": {str(row["dose_mg"]): row["growth_auc"]
                                   for row in ladder_rows},
            "inhibition_by_dose": {str(row["dose_mg"]): row["inhibition_max"]
                                   for row in ladder_rows},
            "cmax_by_dose": {str(row["dose_mg"]): row["cmax_mg_per_l"]
                             for row in ladder_rows},
            "c_max_mg_per_l": float(C_p.max()),
            "conc_auc_mg_h_per_l": float(np.trapezoid(C_p, t_hours)),
            "formal_lemmas_verified": formal["verified"],
        },
    }
    return payload, formal


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=DRIFT_ROOT / "output" / "integration.json")
    ap.add_argument("--provenance", type=Path,
                    default=DRIFT_ROOT / "output" / "regulatory_provenance.json")
    ap.add_argument("--dose", type=float, default=TRIAL_DOSE_MG,
                    help="the trial dose in mg (default: %(default)s)")
    args = ap.parse_args(argv)

    # The record a previous run left: it is the fallback for a
    # repository that is not a git checkout here.  It is read BEFORE the
    # run, and the artefacts are written only after every assertion has
    # passed.
    carried = load_carried_provenance(args.provenance)

    try:
        payload, formal = run_integration(args)
        provenance = build_provenance(payload, formal, carried)
    except IntegrationFailure as failure:
        print(f"INTEGRATION FAILED: {failure}", file=sys.stderr)
        return 1

    # The chain must name four real heads.
    unresolved = [name for name, row in provenance["repos"].items()
                  if row["sha"] == "unknown"]
    if unresolved:
        print(f"INTEGRATION FAILED: the provenance chain does not resolve "
              f"the heads of {unresolved}", file=sys.stderr)
        return 1

    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2) + "\n",
                                encoding="utf-8")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    growths = [row["mean_growth"] for row in payload["dose_ladder"]]
    print(f"INTEGRATION PASSED: C_max = {payload['cmax']:.4f} mg/L, "
          f"max inhibition = {payload['inhibition_max']:.4f}, "
          f"growth monotone in dose ({growths}), "
          f"cell death reached by an infeasible metabolic problem, "
          f"{payload['formal_gate']['verified']} lemmas verified by QED.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
