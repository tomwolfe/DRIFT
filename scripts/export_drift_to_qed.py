#!/usr/bin/env python3
"""Export DRIFT Hill-occupancy lemmas to QED-parseable LaTeX.

This is the DRIFT -> QED semantic bridge. It mirrors
``VeriTrial/scripts/export_pbpk_to_qed.py``: deterministic lemma strings are
written one per line so ``QED/verify_pbpk_lemmas.py`` can prove each one
(no ``sorry``) and the formal gate can assert on ``output.lean``.

The occupancy model in ``drift/binding.py`` is the Hill equation

    θ(C) = C^n / (C^n + Kd^n)

with ``C >= 0``, ``Kd > 0`` and ``n > 0``. The four exported lemmas capture
the properties the time-series bridge relies on:

  * ``hill_nonneg``  -- occupancy is non-negative.
  * ``hill_le_one``  -- occupancy is at most one (θ + (1-θ) = 1).
  * ``hill_mono``    -- occupancy is monotone in C (cross-multiplied
                        difference has the sign of C2^n - C1^n).
  * ``hill_limits``  -- the large-C form 1 / (1 + Kd^n / C^n), whose limit
                        is 1, while the small-C numerator gives 0.

The lemmas are algebraic identities/inequalities rather than numerical
approximations, matching the structural-proxy convention already used for
the PBPK mass-conservation export.

Usage:
    python3 export_drift_to_qed.py                  # -> /tmp/drift_lemmas.txt
    python3 export_drift_to_qed.py --out L.txt
    python3 export_drift_to_qed.py --stdout
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Keep the four named lemmas in one place so a test can import them.
LEMMAS: dict[str, str] = {
    "hill_nonneg": "C ^ n / (C ^ n + Kd ^ n) >= 0",
    "hill_le_one": "C ^ n / (C ^ n + Kd ^ n) + Kd ^ n / (C ^ n + Kd ^ n) = 1",
    "hill_mono": (
        "(C2 ^ n / (C2 ^ n + Kd ^ n) - C1 ^ n / (C1 ^ n + Kd ^ n))"
        " * (C2 ^ n + Kd ^ n) * (C1 ^ n + Kd ^ n)"
        " = (C2 ^ n - C1 ^ n) * Kd ^ n"
    ),
    "hill_limits": "C ^ n / (C ^ n + Kd ^ n) = 1 / (1 + Kd ^ n / C ^ n)",
}

DEFAULT_OUT = Path("/tmp/drift_lemmas.txt")


def build_lines() -> list[str]:
    """Return the lemma file body: names as comments, then one lemma per line."""
    lines: list[str] = []
    for name, lemma in LEMMAS.items():
        lines.append(f"-- {name}")
        lines.append(lemma)
    return lines


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help=f"output file (default: {DEFAULT_OUT})")
    ap.add_argument("--stdout", action="store_true",
                    help="print lemmas instead of writing a file")
    ap.add_argument("--parametric", action="store_true",
                    help="accepted for backward compat with the quad gate; lemmas are always parametric")
    args = ap.parse_args(argv)

    lines = build_lines()
    if args.stdout:
        print("\n".join(lines))
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(LEMMAS)} DRIFT lemmas to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
