"""Bake an ACEScc -> ACEScc negative-film LMT (.cube).

    python scripts/aces_lmt/bake.py kodak_vision3_250d --size 65 --out experiment_results/luts
"""
import argparse, time
from pathlib import Path
import numpy as np
from spektrafilm_lut_creator.aces_lmt import NegativeADXModel, NegativeADXSpec, bake_acescc_lmt, write_cube, write_lut_license_files

ap = argparse.ArgumentParser()
ap.add_argument("films", nargs="+")
ap.add_argument("--size", type=int, default=65)
ap.add_argument("--metric", default="apd", choices=["apd", "print"])
ap.add_argument("--balance", default="grey", choices=["grey", "printer", "aces", "none"])
ap.add_argument("--calibration", default="none", choices=["none", "gamma"])
ap.add_argument("--exposure-ev", type=float, default=0.0)
ap.add_argument("--domain", type=float, nargs=2, default=(0.0, 1.0),
                help="ACEScc input domain; alcedo_studio hard-wires [0, 1]")
ap.add_argument("--out", default="experiment_results/luts")
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
write_lut_license_files(out)
for film in a.films:
    t = time.perf_counter()
    spec = NegativeADXSpec(film_profile=film, density_metric=a.metric, balance=a.balance, calibration=a.calibration, exposure_ev=a.exposure_ev)
    m = NegativeADXModel(spec)
    table = bake_acescc_lmt(m, a.size, domain=tuple(a.domain))
    name = f"spektrafilm_{film}_adx_acescc_lmt_{a.metric}_{a.balance}{'_cal' if a.calibration != 'none' else ''}_{a.size}.cube"
    write_cube(table, out / name, title=f"spektrafilm {film} negative -> ADX16 -> ACEScc LMT", domain=tuple(a.domain), comments=[
        "ACEScc (AP1) -> ACEScc (AP1) look modification transform.",
        f"Scene -> {film} negative (spektrafilm) -> printing density ({a.metric}) -> ADX16 (ST 2065-3)",
        "-> Academy ADX16 IDT -> ACES2065-1 -> ACEScc. No print stock, no display rendering: apply before the DRT.",
        f"balance={a.balance} (18% grey preserved: {a.balance != 'none'}), calibration={a.calibration}, exposure_ev={a.exposure_ev}",
        f"camera CC filtration (log10 E): {np.array2string(m.log_e_offset, precision=4)}",
        f"printer-light CDD offsets: {np.array2string(m.cdd_offset, precision=4)}, ACES gain after IDT: {np.array2string(m.aces_gain, precision=4)}",
        f"per-stock CDD gains (calibration): {np.array2string(m.cdd_gain, precision=4)}",
    ])
    print(f"{name}: {time.perf_counter() - t:.1f}s  out range {table.min():.3f}..{table.max():.3f}")
