"""Compare ACEScc LMT cubes numerically (scene side and after the ACES 2.0 DRT).

    python scripts/aces_lmt/analyze_luts.py --luts a.cube b.cube ... --out experiment_results/analysis
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import colour
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from evaluate import read_cube, apply_trilinear, Drt  # noqa: E402
from spektrafilm_lut_creator.aces_lmt import acescc_encode, acescc_decode, _AP0_TO_AP1, _AP1_TO_AP0  # noqa: E402

AP1 = colour.RGB_COLOURSPACES["ACEScg"]
Y_AP1 = AP1.matrix_RGB_to_XYZ[1]


def lmt(table, ap1):
    if table is None:
        return ap1
    return acescc_decode(apply_trilinear(table, acescc_encode(ap1)))


def to_display(drt, ap1):
    return drt(np.asarray(ap1 @ _AP1_TO_AP0.T))


def display_lab(disp):
    return colour.XYZ_to_Lab(colour.sRGB_to_XYZ(disp), colour.CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["D65"])


def lch(lab):
    return lab[..., 0], np.hypot(lab[..., 1], lab[..., 2]), np.degrees(np.arctan2(lab[..., 2], lab[..., 1])) % 360


def colorchecker_ap1():
    cc = colour.CCS_COLOURCHECKERS["ColorChecker24 - After November 2014"]
    names = list(cc.data.keys())
    xyz = colour.xyY_to_XYZ(np.array(list(cc.data.values())))
    # D50 reference -> ACES white (D60) via CAT02, as an IDT would for a D50-lit chart
    ap1 = colour.XYZ_to_RGB(xyz, AP1, illuminant=cc.illuminant, chromatic_adaptation_transform="CAT02")
    ap1 *= 0.18 / (ap1[21] @ Y_AP1 / 1.0) * (0.19 / 0.19)  # Neutral 5 (~19 %) -> 0.18
    return names, ap1


def oklab_scene(ap1):
    xyz = ap1 @ AP1.matrix_RGB_to_XYZ.T
    return colour.XYZ_to_Oklab(xyz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--luts", nargs="+", required=True)
    ap.add_argument("--out", default="experiment_results/analysis")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    drt = Drt()
    luts = [("ACES 2.0 only", None)] + [(Path(p).stem[:28], read_cube(Path(p))) for p in a.luts]

    stops = np.linspace(-8, 10, 181)
    grey = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
    names, cc = colorchecker_ap1()
    focus = {"dark skin": 0, "light skin": 1, "blue sky": 2, "foliage": 3, "orange": 6, "blue": 12, "green": 13, "red": 14, "yellow": 15, "cyan": 17}
    ref_lab = display_lab(to_display(drt, cc))
    ref_L, ref_C, ref_h = lch(ref_lab)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
    print("\n== neutral axis: scene side (output stops re 0.18) and after DRT (display L*) ==")
    print(f"{'LUT':30s} {'grey->':>7s} {'slope@0':>7s} {'-6':>6s} {'-4':>6s} {'-2':>6s} {'+2':>6s} {'+4':>6s} {'+6':>6s} {'+8':>6s} | L*(grey) {'dispL* -4/+2/+4/+6':>22s} | sys slope@grey")
    for name, t in luts:
        o = lmt(t, grey)
        so = np.log2(np.fmax(o @ Y_AP1, 1e-9) / 0.18)
        disp = to_display(drt, o)
        L = display_lab(disp)[:, 0]
        i = lambda s: np.argmin(abs(stops - s))
        slope = (so[i(0.5)] - so[i(-0.5)])
        # system gamma: dlog(display Y)/dlog(scene) at grey
        Yd = colour.sRGB_to_XYZ(disp)[:, 1]
        sys_slope = (np.log2(Yd[i(0.5)]) - np.log2(Yd[i(-0.5)]))
        print(f"{name:30s} {so[i(0)]:7.2f} {slope:7.2f} " + " ".join(f"{so[i(s)]:6.2f}" for s in (-6, -4, -2, 2, 4, 6, 8))
              + f" | {L[i(0)]:6.1f}   " + "/".join(f"{L[i(s)]:.0f}" for s in (-4, 2, 4, 6)) + f"       | {sys_slope:.2f}")
        axs[0].plot(stops, so, label=name)
        axs[1].plot(stops, L, label=name)
    axs[0].plot(stops, stops, "k:", lw=.8)
    axs[0].set(title="LMT neutral axis (scene side)", xlabel="scene stops re 0.18", ylabel="LMT out stops re 0.18 (Y)")
    axs[1].set(title="System tone curve: LMT + ACES 2.0 SDR", xlabel="scene stops re 0.18", ylabel="display L*")

    print("\n== ColorChecker after DRT (vs ACES 2.0 alone): dL*, C* ratio, dh (deg) ==")
    hdr = "".join(f"{k:>16s}" for k in focus)
    print(f"{'LUT':30s} {'meanCratio':>10s} {'mean|dh|':>8s}" + hdr)
    for name, t in luts[1:]:
        lab = display_lab(to_display(drt, lmt(t, cc)))
        Lx, Cx, hx = lch(lab)
        cr = Cx / np.fmax(ref_C, 1e-6)
        dh = (hx - ref_h + 180) % 360 - 180
        chrom = ref_C > 10
        row = "".join(f"  {Lx[j]-ref_L[j]:+4.0f}/{cr[j]:4.2f}/{dh[j]:+4.0f}" for j in focus.values())
        print(f"{name:30s} {cr[chrom].mean():10.2f} {np.abs(dh[chrom]).mean():8.1f}" + row)

    # saturation vs exposure, scene side, measured in Oklab chroma relative to input
    print("\n== scene-side chroma gain (Oklab C_out/C_in) vs exposure, mean over 12 hues at C=0.08 ==")
    hues = np.radians(np.arange(0, 360, 30))
    ev = np.array([-4, -2, 0, 2, 4, 6])
    for name, t in luts[1:]:
        gains = []
        hue_gain = []
        for e in ev:
            grey_lab = oklab_scene(np.array([[0.18 * 2.0 ** e] * 3]))[0]
            lab = np.stack([np.full(12, grey_lab[0]), 0.08 * np.cos(hues) * (grey_lab[0] / 0.56), 0.08 * np.sin(hues) * (grey_lab[0] / 0.56)], -1)
            xyz = colour.Oklab_to_XYZ(lab)
            ap1 = xyz @ AP1.matrix_XYZ_to_RGB.T
            o = oklab_scene(lmt(t, ap1))
            cin = np.hypot(lab[:, 1], lab[:, 2]); cout = np.hypot(o[:, 1], o[:, 2])
            gains.append((cout / cin).mean())
            if e == 0:
                hue_gain = cout / cin
                hue_shift = (np.degrees(np.arctan2(o[:, 2], o[:, 1]) - hues) + 180) % 360 - 180
        print(f"{name:30s} " + " ".join(f"{e:+d}EV:{g:4.2f}" for e, g in zip(ev, gains))
              + " | @grey per hue(0..330 by 30) gain " + " ".join(f"{g:.2f}" for g in hue_gain)
              + " | dh " + " ".join(f"{d:+.0f}" for d in hue_shift))
        axs[2].plot(ev, gains, marker="o", label=name)
    axs[2].axhline(1, color="k", lw=.8, ls=":")
    axs[2].set(title="Scene-side chroma gain vs exposure", xlabel="EV re grey", ylabel="Oklab C out / C in")
    for ax in axs:
        ax.grid(alpha=.3)
    axs[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "lut_analysis.png", dpi=110)
    print("\nplot:", out / "lut_analysis.png")


if __name__ == "__main__":
    main()
