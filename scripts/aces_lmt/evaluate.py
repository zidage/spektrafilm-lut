"""Render test raws through ACES 2.0 with and without ACEScc LMT cubes.

Emulates alcedo_studio's grade stage: scene-linear AP1 -> ACEScc -> clamp to
[0, 1] -> trilinear 3D LUT -> ACEScc decode -> DRT (here: OCIO ACES 2.0
SDR 100 nit Rec.709 on an sRGB display, studio-config v4.0.0).

    python scripts/aces_lmt/evaluate.py --luts a.cube b.cube --raws r1.CR3 r2.RAF --out experiment_results/renders
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import colour
import numpy as np
import PyOpenColorIO as ocio
import rawpy
from PIL import Image, ImageDraw

from spektrafilm_lut_creator.aces_lmt import acescc_decode, acescc_encode, _AP0_TO_AP1, _AP1_TO_AP0

CONFIG = "studio-config-v4.0.0_aces-v2.0_ocio-v2.5"
# Alcedo Studio default: Rec.709 primaries, gamma 2.2 encoding, 100 nit.
DISPLAY, VIEW = "Gamma 2.2 Rec.709 - Display", "ACES 2.0 - SDR 100 nits (Rec.709)"


def read_cube(path: Path) -> np.ndarray:
    size, rows = None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("LUT_3D_SIZE"):
            size = int(line.split()[1])
        elif re.match(r"^[-+.\d]", line):
            rows.append([float(v) for v in line.split()[:3]])
    table = np.asarray(rows, dtype=np.float64)
    return table.reshape(size, size, size, 3)  # [b, g, r]


def apply_trilinear(table: np.ndarray, cc: np.ndarray) -> np.ndarray:
    """alcedo SampleLut3d: clamp to [0, 1], lattice = u * (N - 1), trilinear."""
    n = table.shape[0]
    u = np.clip(cc, 0.0, 1.0) * (n - 1)
    i0 = np.clip(np.floor(u).astype(np.int64), 0, n - 2)
    f = u - i0
    r0, g0, b0 = i0[..., 0], i0[..., 1], i0[..., 2]
    fr, fg, fb = f[..., 0:1], f[..., 1:2], f[..., 2:3]
    out = 0.0
    for db, wb in ((0, 1 - fb), (1, fb)):
        for dg, wg in ((0, 1 - fg), (1, fg)):
            for dr, wr in ((0, 1 - fr), (1, fr)):
                out = out + table[b0 + db, g0 + dg, r0 + dr] * (wb * wg * wr)
    return out


def load_raw(path: Path, max_side: int) -> np.ndarray:
    with rawpy.imread(str(path)) as raw:
        rgb = raw.postprocess(
            output_color=rawpy.ColorSpace.ACES, gamma=(1, 1), no_auto_bright=True,
            use_camera_wb=True, output_bps=16, half_size=True,
            highlight_mode=rawpy.HighlightMode.Clip,
        )
    img = rgb.astype(np.float32) / 65535.0  # sensor clip -> 1.0 (alcedo convention)
    h, w = img.shape[:2]
    step = max(1, int(np.ceil(max(h, w) / max_side)))
    return img[::step, ::step]


def exposure_to_grey(aces: np.ndarray, target: float = 0.18) -> float:
    """EV that puts the scene's log-average luminance at ``target``."""
    y = aces @ np.array([0.3439664498, 0.7281660966, -0.0721325464])  # AP0 -> Y
    y = np.clip(y, 1e-5, None)
    return float(np.log2(target / np.exp(np.mean(np.log(y)))))


class Drt:
    def __init__(self):
        cfg = ocio.Config.CreateFromBuiltinConfig(CONFIG)
        t = ocio.DisplayViewTransform(src="ACES2065-1", display=DISPLAY, view=VIEW)
        self.cpu = cfg.getProcessor(t).getDefaultCPUProcessor()

    def __call__(self, aces: np.ndarray) -> np.ndarray:
        buf = np.ascontiguousarray(aces.astype(np.float32).reshape(-1, 3))
        self.cpu.applyRGB(buf)
        code = np.clip(buf.reshape(aces.shape).astype(np.float64), 0.0, 1.0)
        # same display light, re-encoded as sRGB for files and colour metrics
        return colour.cctf_encoding(code ** 2.2, function="sRGB")


def with_lmt(aces: np.ndarray, table: np.ndarray) -> np.ndarray:
    cc = acescc_encode(aces @ _AP0_TO_AP1.T)
    cc = apply_trilinear(table, cc)
    return acescc_decode(cc) @ _AP1_TO_AP0.T


def label(img: np.ndarray, text: str) -> Image.Image:
    im = Image.fromarray((img * 255 + 0.5).astype(np.uint8))
    d = ImageDraw.Draw(im)
    d.rectangle([0, 0, 8 + 7 * len(text), 18], fill=(0, 0, 0))
    d.text((4, 3), text, fill=(255, 255, 255))
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--luts", nargs="+", required=True)
    ap.add_argument("--raws", nargs="+", required=True)
    ap.add_argument("--out", default="experiment_results/renders")
    ap.add_argument("--max-side", type=int, default=1200)
    ap.add_argument("--auto-exposure", action="store_true",
                    help="scale each raw so its log-average luminance sits at 0.18")
    ap.add_argument("--ev", type=float, default=0.0, help="extra exposure (EV) for all raws")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    drt = Drt()
    tables = [(Path(p).stem, read_cube(Path(p))) for p in a.luts]
    for rp in a.raws:
        rp = Path(rp)
        try:
            aces = load_raw(rp, a.max_side).astype(np.float64)
        except rawpy.LibRawError as e:
            print(f"skip {rp.name}: {e}")
            continue
        ev = a.ev + (exposure_to_grey(aces) if a.auto_exposure else 0.0)
        aces *= 2.0 ** ev
        panels = [label(drt(aces), f"{rp.name}  no LMT  ({ev:+.2f} EV)")]
        for name, table in tables:
            panels.append(label(drt(with_lmt(aces, table)), name))
        w, h = panels[0].size
        sheet = Image.new("RGB", (w * len(panels), h))
        for k, p in enumerate(panels):
            sheet.paste(p, (k * w, 0))
        dst = out / f"{rp.stem}.jpg"
        sheet.save(dst, quality=90)
        print(dst)


if __name__ == "__main__":
    main()
