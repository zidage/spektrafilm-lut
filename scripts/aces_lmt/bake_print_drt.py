"""Experimental: bake DRT^-1 o spektrafilm(film -> print -> scan) as an ACEScc LMT."""
import argparse, time
from pathlib import Path
from spektrafilm_lut_creator.aces_lmt import PrintDRTModel, PrintDRTSpec, bake_acescc_lmt, write_cube

ap = argparse.ArgumentParser()
ap.add_argument("films", nargs="+")
ap.add_argument("--print", dest="print_profile", default=None)
ap.add_argument("--size", type=int, default=65)
ap.add_argument("--no-white", action="store_true")
ap.add_argument("--out", default="experiment_results/luts_print")
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
for film in a.films:
    t = time.perf_counter()
    m = PrintDRTModel(PrintDRTSpec(film_profile=film, print_profile=a.print_profile,
                                   white_correction=not a.no_white, black_correction=not a.no_white))
    table = bake_acescc_lmt(m, a.size, domain=(0.0, 1.0))
    pr = m._pipeline.print.info.stock
    name = f"spektrafilm_{film}_{pr}{'_nowhite' if a.no_white else ''}_invACES2_acescc_lmt_{a.size}.cube"
    write_cube(table, out / name, title=f"spektrafilm {film} > {pr} (inverse ACES 2.0 SDR) ACEScc LMT", domain=(0.0, 1.0),
               comments=["EXPERIMENTAL: DRT^-1 o spektrafilm film->print->scan. Use only with ACES 2.0 SDR 100 nit Rec.709."])
    print(f"{name}: {time.perf_counter()-t:.1f}s  range {table.min():.3f}..{table.max():.3f}")
