"""Measure spektrafilm's own display rendering (film -> print -> scan, as in the GUI)
with the same metrics as analyze_luts.py (neutral tone curve, ColorChecker L*C*h)."""
from __future__ import annotations

import sys
from pathlib import Path

import colour
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from analyze_luts import colorchecker_ap1, display_lab, lch, to_display, Y_AP1  # noqa: E402
from evaluate import Drt  # noqa: E402
from spektrafilm_lut_creator.aces_lmt import _AP1_TO_AP0  # noqa: E402


def spektrafilm_display(film, print_stock):
    from spektrafilm.runtime.params_builder import digest_params, init_params
    from spektrafilm.runtime.pipeline import SimulationPipeline
    p = init_params(film_profile=film, print_profile=print_stock)
    p.debug.lut_mode = True
    p.io.input_color_space = "ACES2065-1"
    p.io.input_cctf_decoding = False
    p.io.output_color_space = "sRGB"
    p.io.output_cctf_encoding = True
    pipe = SimulationPipeline(digest_params(p))

    def f(ap1):
        aces = np.asarray(ap1) @ _AP1_TO_AP0.T
        out = pipe.process(aces.reshape(1, -1, 3).astype(np.float64))
        return np.clip(np.asarray(out).reshape(aces.shape), 0, 1)
    return f


def main():
    drt = Drt()
    stops = np.linspace(-8, 10, 181)
    grey = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
    names, cc = colorchecker_ap1()
    focus = {"dark skin": 0, "light skin": 1, "blue sky": 2, "foliage": 3, "orange": 6, "blue": 12, "green": 13, "red": 14, "yellow": 15, "cyan": 17}
    ref_L, ref_C, ref_h = lch(display_lab(to_display(drt, cc)))
    i = lambda s: np.argmin(abs(stops - s))
    chains = [("ACES 2.0 only", lambda x: to_display(drt, x))]
    for film, pr in (("kodak_vision3_250d", "kodak_2383"), ("kodak_portra_400", "kodak_portra_endura"),
                     ("fujifilm_pro_400h", "fujifilm_crystal_archive_typeii")):
        chains.append((f"sf {film[:14]}>{pr[:10]}", spektrafilm_display(film, pr)))
    print(f"{'display rendering':34s} L*grey  L* -6/-4/-2/+2/+4/+6/+8  sys slope@grey | meanCratio mean|dh| " + "".join(f"{k:>16s}" for k in focus))
    for name, f in chains:
        disp = f(grey)
        L = display_lab(disp)[:, 0]
        Yd = colour.sRGB_to_XYZ(disp)[:, 1]
        sl = np.log2(Yd[i(0.5)]) - np.log2(Yd[i(-0.5)])
        Lx, Cx, hx = lch(display_lab(f(cc)))
        cr = Cx / np.fmax(ref_C, 1e-6)
        dh = (hx - ref_h + 180) % 360 - 180
        m = ref_C > 10
        row = "".join(f"  {Lx[j]-ref_L[j]:+4.0f}/{cr[j]:4.2f}/{dh[j]:+4.0f}" for j in focus.values())
        print(f"{name:34s} {L[i(0)]:5.1f}  " + "/".join(f"{L[i(s)]:.0f}" for s in (-6, -4, -2, 2, 4, 6, 8))
              + f"   {sl:5.2f}        | {cr[m].mean():6.2f} {np.abs(dh[m]).mean():6.1f}" + row)


if __name__ == "__main__":
    main()
