"""Experimental: bake DRT^-1 o spektrafilm(film -> print -> scan) as an ACEScc LMT."""
import argparse, time
from pathlib import Path
from spektrafilm_lut_creator.aces_lmt import PrintDRTModel, PrintDRTSpec, bake_acescc_lmt, write_cube, write_lut_license_files

PRINT_SHORT = {"kodak_2383": "2383", "kodak_2393": "2393", "kodak_portra_endura": "Endura",
               "kodak_endura_premier": "EnduraPremier", "kodak_ektacolor_edge": "Edge",
               "kodak_supra_endura": "Supra", "kodak_ultra_endura": "Ultra",
               "fujifilm_crystal_archive_typeii": "CA", "slide": "Slide"}
FILM_SHORT = {"vision3_50d": "Vision3-50D", "vision3_200t": "Vision3-200T", "vision3_250d": "Vision3-250D",
              "vision3_500t": "Vision3-500T", "verita_200d": "Verita200D", "portra_160": "Portra160",
              "portra_400": "Portra400", "portra_800": "Portra800", "portra_800_push1": "Portra800P1",
              "portra_800_push2": "Portra800P2", "ektar_100": "Ektar100", "gold_200": "Gold200",
              "ultramax_400": "UltraMax400", "c200": "C200", "pro_400h": "Pro400H", "xtra_400": "Xtra400",
              "velvia_100": "Velvia100", "provia_100f": "Provia100F", "ektachrome_100": "Ektachrome100",
              "kodachrome_64": "Kodachrome64"}


def short_name(film: str, print_stock: str, neutralize: bool, hue: float, chroma: float) -> str:
    """Brand_Film_Print[_calib], e.g. Kodak_Vision3-250D_2383_NH."""
    brand, _, stock = film.partition("_")
    brand = {"kodak": "Kodak", "fujifilm": "Fuji"}.get(brand, brand.capitalize())
    calib = ("N" if neutralize else "") + ("H" if hue else "") + ("C" if chroma != 1 else "")
    parts = [brand, FILM_SHORT.get(stock, stock), PRINT_SHORT.get(print_stock, print_stock)]
    return "_".join(parts + ([calib] if calib else []))


ap = argparse.ArgumentParser()
ap.add_argument("films", nargs="+")
ap.add_argument("--print", dest="print_profile", default=None)
ap.add_argument("--size", type=int, default=65)
ap.add_argument("--no-white", action="store_true")
ap.add_argument("--neutralize", action="store_true", help="remove the print's neutral-scale crossover")
ap.add_argument("--hue-preserve", type=float, default=0.0, help="0 = film hues, 1 = DRT hues")
ap.add_argument("--chroma", type=float, default=1.0)
ap.add_argument("--long-names", action="store_true", help="old descriptive file names")
ap.add_argument("--out", default="experiment_results/luts_print")
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
write_lut_license_files(out)
for film in a.films:
    t = time.perf_counter()
    m = PrintDRTModel(PrintDRTSpec(film_profile=film, print_profile=a.print_profile,
                                   white_correction=not a.no_white, black_correction=not a.no_white,
                                   neutralize=a.neutralize, hue_preserve=a.hue_preserve, chroma_gain=a.chroma))
    table = bake_acescc_lmt(m, a.size, domain=(0.0, 1.0))
    pr = "slide" if m.reversal else m._pipeline.print.info.stock
    tag = ("_nowhite" if a.no_white else "") + ("_neutral" if a.neutralize else "")         + (f"_hue{a.hue_preserve:g}" if a.hue_preserve else "") + (f"_chroma{a.chroma:g}" if a.chroma != 1 else "")
    name = f"spektrafilm_{film}_{pr}{tag}_invACES2_acescc_lmt_{a.size}.cube"
    if not a.long_names:
        name = short_name(film, pr, a.neutralize, a.hue_preserve, a.chroma) + ".cube"
    write_cube(table, out / name, title=name[:-5] if not a.long_names else f"spektrafilm {film} > {pr} ACEScc LMT (for ACES 2.0 SDR)", domain=(0.0, 1.0),
               comments=[("ACES2.0 SDR^-1 o spektrafilm reversal film -> scan." if m.reversal else "ACES2.0 SDR^-1 o spektrafilm film -> print -> scan.") + " Use only with ACES 2.0 SDR 100 nit, Rec.709 primaries, gamma 2.2 display encoding.",
                         f"neutralize={a.neutralize} hue_preserve={a.hue_preserve} chroma_gain={a.chroma}"])
    print(f"{name}: {time.perf_counter()-t:.1f}s  range {table.min():.3f}..{table.max():.3f}")
