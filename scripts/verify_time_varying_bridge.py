#!/usr/bin/env python3
"""Integration check: VeriTrial-style C(t) -> DRIFT time-varying simulation.

Builds a genuinely time-varying concentration trace (``1 + 0.5*sin(t)``),
feeds it through :meth:`drift.workbench.Workbench.from_veritrial`, runs a
small Monte Carlo ensemble and asserts that drug exposure suppresses growth
below the basal rate.  Also checks that the audit log records the array
(interpolated) concentration source, so a scalar placeholder cannot pass.

Exits non-zero on any failure.  Intended as a tether verification command.
"""

from __future__ import annotations

import numpy as np

from drift.workbench import Workbench


def main() -> int:
    t_eval = np.linspace(0.0, 24.0, 49)
    c_p = 1.0 + 0.5 * np.sin(np.linspace(0.0, 2.0 * np.pi, 49))
    assert np.ptp(c_p) > 0.5, "test trace must actually vary"

    wb = Workbench.from_veritrial(
        {"t_eval": t_eval, "C_p": c_p, "patient_id": "integration"},
        steps=60,
        model_name="textbook",
    )
    if not hasattr(wb, "run_monte_carlo"):
        raise SystemExit("from_veritrial did not return a configured Workbench")

    result = wb.run_monte_carlo(n_sims=2, steps=60, n_jobs=1)
    basal = float(result["basal_growth"])
    growths = [float(np.mean(h["growth"])) for h in result["histories"]]
    sources = [h.get("concentration_source") for h in result["histories"]]
    if not growths:
        raise SystemExit("Monte Carlo returned no histories")
    if not all(g < basal for g in growths):
        raise SystemExit(
            f"time-varying exposure did not suppress growth: basal={basal} "
            f"growths={growths}")
    if any(src != "array_interp" for src in sources):
        raise SystemExit(
            f"expected array_interp concentration source, got {sources}")
    print(
        f"time-varying bridge OK: basal={basal:.6f} "
        f"growths={[round(g, 6) for g in growths]} sources={sources}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
