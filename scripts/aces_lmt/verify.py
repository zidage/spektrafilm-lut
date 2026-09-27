"""Quantitative checks for an ACEScc LMT cube.

1. Tone scale: neutral ramp -40..+10 stops through the cube (alcedo sampling)
   vs the exact model; mid-grey error, per-channel crossover, output range.
2. LUT fidelity on real raws: trilinear cube vs exact model, measured after
   the ACES 2.0 DRT as CIE dE2000 on display sRGB.
"""
import argparse, sys
from pathlib import Path
import numpy as np, colour
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from evaluate import read_cube, apply_trilinear, load_raw, exposure_to_grey, Drt, with_lmt
from spektrafilm_lut_creator.aces_lmt import NegativeADXModel, NegativeADXSpec, acescc_encode, acescc_decode, _AP0_TO_AP1, _AP1_TO_AP0

ap = argparse.ArgumentParser()
ap.add_argument("cube"); ap.add_argument("--film", required=True)
ap.add_argument("--metric", default="apd"); ap.add_argument("--calibration", default="none"); ap.add_argument("--raws", nargs="*", default=[])
ap.add_argument("--out", default="experiment_results/plots")
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
table = read_cube(Path(a.cube))
model = NegativeADXModel(NegativeADXSpec(film_profile=a.film, density_metric=a.metric, calibration=a.calibration))

stops = np.linspace(-12, 12, 97)
aces = (0.18 * 2 ** stops)[:, None] * np.ones(3)
lut = with_lmt(aces, table)
exact = model.aces_out(aces)
ap1_lut = lut @ _AP0_TO_AP1.T
s_lut = np.log2(np.fmax(ap1_lut, 1e-9) / 0.18)
s_exact = np.log2(np.fmax(exact @ _AP0_TO_AP1.T, 1e-9) / 0.18)
i0 = np.argmin(abs(stops))
print(f"mid-grey out (AP1 stops re 0.18): lut {s_lut[i0].round(4)}  exact {s_exact[i0].round(4)}")
for s in (-8, -6, -4, -2, 2, 4, 6, 8, 10):
    i = np.argmin(abs(stops - s)); print(f"  scene {s:+3d} st -> out {s_lut[i].round(2)}  crossover(max-min) {np.ptp(s_lut[i]):.3f} st")
print(f"output linear max {ap1_lut.max():.2f} ({np.log2(ap1_lut.max()/0.18):+.2f} st), input domain covers "
      f"{np.log2(acescc_decode(0.0)/0.18):+.1f}..{np.log2(acescc_decode(1.0)/0.18):+.1f} st")
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(stops, stops, "k:", lw=1, label="identity")
for c, col in enumerate("rgb"):
    ax.plot(stops, s_exact[:, c], color=col, lw=1, alpha=.5)
    ax.plot(stops, s_lut[:, c], color=col, ls="--", lw=1.2, label=f"{col.upper()} (65³ cube)" if c == 0 else None)
ax.axvspan(np.log2(acescc_decode(0.0)/0.18), np.log2(acescc_decode(1.0)/0.18), color="0.9", zorder=-1, label="cube domain ACEScc [0,1]")
ax.set(xlabel="scene exposure (stops re 18% grey)", ylabel="LMT output (stops re 0.18, AP1)", title=f"{a.film} → ADX16 → ACES LMT tone scale")
ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(out / f"tonescale_{a.film}.png", dpi=120)

if a.raws:
    drt = Drt(); des = []
    for rp in a.raws:
        try: img = load_raw(Path(rp), 600).astype(np.float64)
        except Exception as e: print("skip", rp, e); continue
        img *= 2 ** exposure_to_grey(img)
        d_lut = drt(with_lmt(img, table))
        d_ex = drt(model.aces_out(img.reshape(-1, 3)).reshape(img.shape))
        lab = lambda x: colour.XYZ_to_Lab(colour.sRGB_to_XYZ(x))
        de = colour.delta_E(lab(d_lut), lab(d_ex), method="CIE 2000")
        des.append(de.ravel()); print(f"  {Path(rp).name}: dE2000 mean {de.mean():.3f}  p99 {np.percentile(de,99):.3f}  max {de.max():.2f}")
    de = np.concatenate(des); print(f"ALL: dE2000 mean {de.mean():.3f} p99 {np.percentile(de,99):.3f} max {de.max():.2f}")
