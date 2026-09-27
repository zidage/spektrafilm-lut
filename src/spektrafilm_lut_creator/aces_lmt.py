"""ACEScc -> ACEScc film-negative LMT (Look Modification Transform).

Models the motion-picture scanning workflow instead of a display-referred
film + print + scan chain:

    scene (ACES2065-1, treated as the real world)
      -> camera exposure onto a virtual colour negative   (spektrafilm filming stage)
      -> development: CMY dye densities                   (spektrafilm develop, ``cmy_film`` tap)
      -> spectral transmittance of the processed negative (channel dye spectra + base / mask)
      -> Academy Printing Density (APD, SMPTE ST 2065-2)  (spectral integration)
      -> ADX16 encoding (SMPTE ST 2065-3)                 (Dmin-subtracted, channel gains)
      -> Academy ADX16 IDT                                 (OCIO builtin ADX16_to_ACES2065-1)
      -> ACES2065-1 -> ACEScc

There is no print stock and no display rendering here: the negative is
scanned to printing density (scene-referred, like a real film scan) and
re-enters ACES, so the full scene dynamic range the negative can record
survives and the DRT downstream (e.g. ACES 2.0) does the tone mapping.

The APD spectral responsivities are SMPTE data; they are downloaded from
pub.smpte.org on first use and cached (not redistributed in this repo).
"""
from __future__ import annotations

import io
import os
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# SMPTE ST 2065-3 ADX16 encoding: CV = gain_c * (APD_c - APD_c,Dmin) * 8000 + 1520.
ADX16_SCALE = 8000.0
ADX16_OFFSET = 1520.0
ADX_CHANNEL_GAINS = np.array([1.00, 0.92, 0.95])
# Laboratory aim density: Cineon LAD 445 -> (445 - 95) / 500 = 0.70 channel
# dependent density above Dmin; the Academy ADX IDT maps a neutral CDD of 0.70
# to ACES 0.18 (REF_PT in ACEScsc.Academy.ADX16_to_ACES.ctl).
LAD_CDD = 0.70

# ACEScc (S-2014-003) encoding bounds.
ACESCC_MIN = (np.log2(2.0 ** -16) + 9.72) / 17.52  # value for linear <= 0
ACESCC_MAX = (np.log2(65504.0) + 9.72) / 17.52     # half-float max

_APD_ZIP_URL = "https://pub.smpte.org/pub/st2065-2/st2065-2-2020.zip"
_APD_CSV_NAME = "st2065-2a-2020.csv"
_APD_INFLUX_CSV_NAME = "st2065-2b-2020.csv"


# ---------------------------------------------------------------------------
# APD spectral responsivities (SMPTE ST 2065-2)


def _apd_cache_dir() -> Path:
    env = os.environ.get("SPEKTRAFILM_APD_DIR")
    if env:
        return Path(env)
    return Path.home() / ".cache" / "spektrafilm" / "smpte_st2065_2"


def _smpte_csv(name: str, cache_dir: Path | None) -> np.ndarray:
    cache_dir = Path(cache_dir) if cache_dir is not None else _apd_cache_dir()
    csv_path = cache_dir / name
    if not csv_path.exists():
        cache_dir.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(_APD_ZIP_URL, timeout=60) as resp:
            payload = resp.read()
        with zipfile.ZipFile(io.BytesIO(payload)) as zf:
            for member in (_APD_CSV_NAME, _APD_INFLUX_CSV_NAME):
                (cache_dir / member).write_bytes(zf.read(member))
    data = np.genfromtxt(csv_path, delimiter=",", skip_header=6)
    return data[~np.isnan(data[:, 0])]


def load_apd_responsivities(cache_dir: Path | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(wavelengths_nm, Pi_APD[n, 3])`` at the standard's 2 nm grid.

    Downloads the ST 2065-2 data supplement from pub.smpte.org the first time
    and caches it under ``cache_dir`` (or ``$SPEKTRAFILM_APD_DIR``).
    """
    data = _smpte_csv(_APD_CSV_NAME, cache_dir)
    return data[:, 0], data[:, 1:4]


def load_apd_influx(cache_dir: Path | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(wavelengths_nm, S_APD[n])``: the ST 2065-2 Annex A printer
    influx (Bell & Howell Model C lamp house + Wratten 2B), 1.0 at 560 nm."""
    data = _smpte_csv(_APD_INFLUX_CSV_NAME, cache_dir)
    return data[:, 0], data[:, 1]


def printing_density(density_spectral: np.ndarray, wavelengths: np.ndarray,
                     apd_wavelengths: np.ndarray, apd_responsivities: np.ndarray) -> np.ndarray:
    """ST 2065-2 eq. (1): APD_c = -log10( sum(Pi_c * T) / sum(Pi_c) ).

    ``density_spectral`` is ``(..., n_wl)`` optical density on ``wavelengths``;
    it is interpolated onto the APD 2 nm grid (the standard requires
    interpolation when sampling differs).
    """
    flat = density_spectral.reshape(-1, density_spectral.shape[-1])
    # np.interp holds the edge value outside the film data range (380-780 nm);
    # the APD responsivities are ~0 below 380 nm.
    resampled = np.empty((flat.shape[0], apd_wavelengths.size))
    for k, wl in enumerate(apd_wavelengths):
        j = np.searchsorted(wavelengths, wl)
        if j <= 0:
            resampled[:, k] = flat[:, 0]
        elif j >= wavelengths.size:
            resampled[:, k] = flat[:, -1]
        else:
            t = (wl - wavelengths[j - 1]) / (wavelengths[j] - wavelengths[j - 1])
            resampled[:, k] = (1 - t) * flat[:, j - 1] + t * flat[:, j]
    transmittance = 10.0 ** (-np.nan_to_num(resampled, nan=0.0))
    weights = apd_responsivities / apd_responsivities.sum(axis=0, keepdims=True)
    apd = -np.log10(np.fmax(transmittance @ weights, 1e-30))
    return apd.reshape(*density_spectral.shape[:-1], 3)


# ---------------------------------------------------------------------------
# OCIO helpers (ACES standard transforms, used as the reference implementation)


def _ocio_cpu(style: str, inverse: bool = False):
    import PyOpenColorIO as ocio
    cfg = ocio.Config.CreateRaw()
    tr = ocio.BuiltinTransform(style)
    if inverse:
        tr.setDirection(ocio.TRANSFORM_DIR_INVERSE)
    return cfg.getProcessor(tr).getDefaultCPUProcessor()


def _apply_cpu(cpu, rgb: np.ndarray) -> np.ndarray:
    buf = np.ascontiguousarray(np.asarray(rgb, dtype=np.float32).reshape(-1, 3))
    cpu.applyRGB(buf)
    return buf.reshape(np.shape(rgb)).astype(np.float64)


# ACEScc (S-2014-003) with the linear extension below zero used by
# alcedo_studio (``cuda_acescc.cuh``): v < 0 encodes as ACESCC_MIN + v so
# small negative AP1 values (e.g. from the IDT's CDD->CID matrix on extreme
# blues) round-trip instead of clamping.
_AP0_TO_AP1 = np.array([
    [1.4514393161, -0.2365107469, -0.2149285693],
    [-0.0765537734, 1.1762296998, -0.0996759264],
    [0.0083161484, -0.0060324498, 0.9977163014],
])
_AP1_TO_AP0 = np.linalg.inv(_AP0_TO_AP1)


def acescc_encode(lin_ap1: np.ndarray) -> np.ndarray:
    v = np.asarray(lin_ap1, dtype=np.float64)
    small = (np.log2(2.0 ** -16 + np.fmax(v, 0.0) * 0.5) + 9.72) / 17.52
    big = (np.log2(np.fmax(v, 2.0 ** -15)) + 9.72) / 17.52
    return np.where(v < 0.0, ACESCC_MIN + v, np.where(v < 2.0 ** -15, small, big))


def acescc_decode(cc: np.ndarray) -> np.ndarray:
    v = np.asarray(cc, dtype=np.float64)
    lin = np.exp2(v * 17.52 - 9.72)
    small = (lin - 2.0 ** -16) * 2.0
    return np.where(v < ACESCC_MIN, v - ACESCC_MIN,
                    np.where(v <= (9.72 - 15.0) / 17.52, small, lin))


def acescc_to_aces(acescc: np.ndarray) -> np.ndarray:
    return acescc_decode(acescc) @ _AP1_TO_AP0.T


def aces_to_acescc(aces: np.ndarray) -> np.ndarray:
    return acescc_encode(np.asarray(aces, dtype=np.float64) @ _AP0_TO_AP1.T)


def adx16_to_aces(adx16_cv: np.ndarray) -> np.ndarray:
    """Academy ADX16 IDT; ``adx16_cv`` in 16-bit code values (float, unclamped)."""
    return _apply_cpu(_ocio_cpu("ADX16_to_ACES2065-1"), np.asarray(adx16_cv) / 65535.0)


def aces_to_adx16(aces: np.ndarray) -> np.ndarray:
    return _apply_cpu(_ocio_cpu("ADX16_to_ACES2065-1", inverse=True), aces) * 65535.0


def _fill_edges(spectra: np.ndarray) -> np.ndarray:
    """Hold the nearest valid value across NaN gaps in a spectral table.

    Profile dye spectra are undefined at the ends (e.g. 380-400 nm for the
    Vision3 stocks) but printing responsivities are not zero there: treating
    NaN as zero density would make that band perfectly transparent and cap
    every printing density with a spurious flare term.
    """
    out = np.array(spectra, dtype=float, copy=True)
    cols = out.reshape(out.shape[0], -1)
    idx = np.arange(cols.shape[0])
    for c in range(cols.shape[1]):
        ok = ~np.isnan(cols[:, c])
        if ok.any() and not ok.all():
            cols[:, c] = np.interp(idx, idx[ok], cols[ok, c])
    return np.nan_to_num(out)


# ---------------------------------------------------------------------------
# The negative -> ADX model


@dataclass
class NegativeADXSpec:
    film_profile: str = "kodak_vision3_250d"
    # Camera exposure relative to the film's rated speed (EV). 0 = box speed.
    exposure_ev: float = 0.0
    # Printing-density metric used by the virtual scanner:
    # 'apd'   -- SMPTE ST 2065-2 Academy Printing Density responsivities (strict).
    # 'print' -- RP 180 printing density of spektrafilm's own printer: the
    #            stock's target print film sensitivity under the neutral-filtered
    #            enlarger lamp.  The stock profiles are fitted against exactly
    #            this printer, so it is the self-consistent metric of the model
    #            (ST 2065-2 notes APD itself is modelled on these print films).
    density_metric: str = "apd"
    # Scanner calibration puts Dmin on the ADX aim (1520) in every mode.  Then
    # an exposed 18 % grey card is balanced to neutral 0.18 by:
    # 'grey'    -- camera CC filtration: per-layer log-exposure offsets so grey
    #              lands on equal CDD in all channels, then one scalar ACES
    #              exposure gain.  Dmin stays on the aim and every layer keeps
    #              its own toe/shoulder, so the (very steep) toe of the Academy
    #              IDT near Dmin is not disturbed.
    # 'printer' -- printer lights: per-channel density offsets to equal CDD,
    #              then a scalar gain.  Moves each channel's Dmin off the aim,
    #              which the IDT toe turns into coloured deep shadows.
    # 'aces'    -- per-channel ACES gain after the IDT (a colourist's offset in
    #           ACEScc) returns grey to 0.18.  Fine for well balanced stocks
    #           (Vision3); pushes saturated colours negative on stocks whose
    #           grey is far from neutral in printing density.
    # 'none'    -- Dmin calibration only; the stock's own grey density and colour
    #           balance go straight through the Academy IDT.
    balance: str = "grey"
    # Per-stock scanner calibration (ST 2065-2 note 3 allows product-specific
    # transforms): 'none', or 'gamma' -- per-channel CDD gains equalising the
    # mid-scale (-2..+2 stop) channel gammas to their mean.  Removes the linear
    # crossover of stocks the Academy IDT was not designed for (still-photo
    # negatives) while keeping overall contrast, toe and shoulder.
    calibration: str = "none"
    # Input gamut compression of the spectral upsampling (runtime default).
    input_gamut_compress: object | None = None
    rgb_to_raw_method: str = "hanatos2025"
    extra: dict = field(default_factory=dict)


class NegativeADXModel:
    """Scene-linear ACES2065-1 -> negative printing density -> ADX16 -> ACES."""

    def __init__(self, spec: NegativeADXSpec):
        from spektrafilm.runtime.params_builder import digest_params, init_params
        from spektrafilm.runtime.pipeline import SimulationPipeline

        from spektrafilm.profiles.io import load_profile

        self.spec = spec
        self.print_stock = load_profile(spec.film_profile).info.target_print or "kodak_2383"
        params = init_params(film_profile=spec.film_profile, print_profile=self.print_stock)
        params.debug.lut_mode = True
        params.io.input_color_space = "ACES2065-1"
        params.io.input_cctf_decoding = False
        params.io.output_cctf_encoding = False
        params.io.scan_film = True
        params.settings.rgb_to_raw_method = spec.rgb_to_raw_method
        if spec.input_gamut_compress is not None:
            params.io.input_gamut_compress = spec.input_gamut_compress
        params = digest_params(params)
        # lut_mode forces camera exposure to 0; the LMT's exposure offset is
        # applied to the scene values instead (identical: raw is linear in E).
        self._pipeline = SimulationPipeline(params)
        film = self._pipeline.film
        self.film = film
        self._wl = np.asarray(film.data.wavelengths, dtype=float)
        self._channel_density = _fill_edges(np.asarray(film.data.channel_density, dtype=float))
        self._base_density = _fill_edges(np.asarray(film.data.base_density, dtype=float))
        self.resp_wl, self.resp = self._responsivities(spec.density_metric, film)

        # Dmin: processed but unexposed negative (cmy_film == 0 -> base + mask).
        self.pd_dmin = self.pd_from_cmy(np.zeros((1, 1, 3)))[0, 0]
        self.log_e_offset = np.zeros(3)  # camera CC filtration (log10 exposure)
        self.cdd_gain = np.ones(3)       # per-stock scanner calibration
        self.cdd_offset = np.zeros(3)    # printer lights
        self.aces_gain = np.ones(3)      # post-IDT gain
        if spec.balance not in ("grey", "printer", "aces", "none"):
            raise ValueError(f"unknown balance {spec.balance!r}")
        if spec.calibration not in ("gamma", "none"):
            raise ValueError(f"unknown calibration {spec.calibration!r}")

        grey = np.full((1, 3), 0.18)
        for _ in range(3 if spec.calibration == "gamma" and spec.balance == "grey" else 1):
            if spec.calibration == "gamma":
                self.cdd_gain = np.ones(3)
                lo, hi = self.cdd(np.array([[0.18 / 4] * 3, [0.18 * 4] * 3]))
                self.cdd_gain = (hi - lo).mean() / (hi - lo)
            if spec.balance == "grey":
                self._solve_cc_filtration(grey)
        if spec.balance == "printer":
            grey_cdd = self.cdd(grey)[0]
            self.cdd_offset = grey_cdd.mean() - grey_cdd
        if spec.balance in ("grey", "printer"):
            self.aces_gain = np.full(3, 0.18 / self.aces_out(grey)[0].mean())
        elif spec.balance == "aces":
            self.aces_gain = 0.18 / self.aces_out(grey)[0]

    def _solve_cc_filtration(self, grey: np.ndarray, iters: int = 8, h: float = 0.02) -> None:
        """Newton-solve per-layer log-exposure offsets so grey has equal CDD."""
        for _ in range(iters):
            c0 = self.cdd(grey)[0]
            err = c0 - c0.mean()
            if np.max(np.abs(err)) < 1e-5:
                break
            base = self.log_e_offset.copy()
            slope = np.empty(3)
            for k in range(3):
                self.log_e_offset = base.copy()
                self.log_e_offset[k] += h
                slope[k] = (self.cdd(grey)[0][k] - c0[k]) / h
            # drive every channel to the mean of the current CDDs
            self.log_e_offset = base - err / np.fmax(slope, 1e-3)

    def _responsivities(self, metric, film):
        apd_wl, apd_pi = load_apd_responsivities()
        if metric == "apd":
            return apd_wl, apd_pi
        if metric == "print":
            from spektrafilm.model.illuminants import standard_illuminant
            pipe = self._pipeline
            sens = np.nan_to_num(10.0 ** np.asarray(pipe.print.data.log_sensitivity, dtype=float))
            lamp = pipe._enlarger_service.enlarger_neutral_illuminant(
                standard_illuminant(pipe.enlarger.illuminant))
            return self._wl, sens * np.asarray(lamp, dtype=float)[:, None]
        raise ValueError(f"unknown density_metric {metric!r}")

    # -- stages -------------------------------------------------------------

    def cmy_film(self, aces: np.ndarray) -> np.ndarray:
        img = np.asarray(aces, dtype=np.float64)
        shape = img.shape
        img = img.reshape(1, -1, 3) * 2.0 ** self.spec.exposure_ev
        log_e = self._pipeline.process(np.fmax(img, 0.0), collect="log_e_film")
        log_e = np.asarray(log_e, dtype=np.float64) + self.log_e_offset
        cmy = self._pipeline.process(log_e, inject="log_e_film", collect="cmy_film")
        return np.asarray(cmy, dtype=np.float64).reshape(shape)

    def pd_from_cmy(self, cmy: np.ndarray) -> np.ndarray:
        """Printing density of the processed negative for dye amounts ``cmy``."""
        dens = np.einsum("...k,lk->...l", cmy, self._channel_density) + self._base_density
        return printing_density(dens, self._wl, self.resp_wl, self.resp)

    def cdd(self, aces: np.ndarray) -> np.ndarray:
        """Channel-dependent density: Dmin subtracted, ADX channel gains."""
        pd = self.pd_from_cmy(self.cmy_film(aces))
        return self.cdd_gain * ADX_CHANNEL_GAINS * (pd - self.pd_dmin)

    def adx16(self, aces: np.ndarray) -> np.ndarray:
        """Float ADX16 code values (unquantised, unclamped)."""
        return (self.cdd(aces) + self.cdd_offset) * ADX16_SCALE + ADX16_OFFSET

    def aces_out(self, aces: np.ndarray) -> np.ndarray:
        return adx16_to_aces(self.adx16(aces)) * self.aces_gain

    def acescc_lmt(self, acescc: np.ndarray) -> np.ndarray:
        return aces_to_acescc(self.aces_out(acescc_to_aces(acescc)))


# ---------------------------------------------------------------------------
# LUT baking


def bake_acescc_lmt(model: NegativeADXModel, size: int = 65,
                    domain: tuple[float, float] = (ACESCC_MIN, ACESCC_MAX),
                    chunk: int = 65 ** 2 * 8) -> np.ndarray:
    """Sample ``model.acescc_lmt`` on a ``size``^3 grid over ``domain``.

    Returns ``(size, size, size, 3)`` indexed ``[b, g, r]`` (Adobe .cube order
    when flattened C-order: R fastest).
    """
    axis = np.linspace(domain[0], domain[1], size)
    b, g, r = np.meshgrid(axis, axis, axis, indexing="ij")
    grid = np.stack([r.ravel(), g.ravel(), b.ravel()], axis=-1)
    out = np.empty_like(grid)
    for s in range(0, grid.shape[0], chunk):
        out[s:s + chunk] = model.acescc_lmt(grid[s:s + chunk])
    return out.reshape(size, size, size, 3)


def write_cube(table: np.ndarray, path: Path, *, title: str,
               domain: tuple[float, float], comments: list[str] = ()) -> None:
    n = table.shape[0]
    lines = [f"# {c}" if c else "#" for c in comments]
    lines.append(f'TITLE "{title}"')
    lines.append("DOMAIN_MIN " + " ".join(f"{domain[0]:.10g}" for _ in range(3)))
    lines.append("DOMAIN_MAX " + " ".join(f"{domain[1]:.10g}" for _ in range(3)))
    lines.append(f"LUT_3D_SIZE {n}")
    for rgb in table.reshape(-1, 3):
        lines.append(" ".join(f"{v:.8f}" for v in rgb))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Experimental: full print chain made scene-referred through an inverse DRT
#
# LMT = DRT^-1 o scan(print(negative(scene))).  Applied before the same DRT the
# system reproduces spektrafilm's display (print) rendering, which is the
# structure of camera "film simulation" LMTs converted in Resolve.  It is tied
# to the DRT it inverts (default: ACES 2.0 SDR 100 nit Rec.709 on sRGB).

_OCIO_CONFIG = "studio-config-v4.0.0_aces-v2.0_ocio-v2.5"


@dataclass
class PrintDRTSpec:
    film_profile: str = "kodak_vision3_250d"
    print_profile: str | None = None  # default: the film's target print
    exposure_ev: float = 0.0
    white_correction: bool = True     # paper white -> display white_level
    black_correction: bool = True     # print Dmax -> display black_level
    white_level: float = 0.98
    black_level: float = 0.005
    # The runtime's output gamut compression also rolls lightness off above
    # 0.7 (cam16ucs lightness_compression); with a DRT downstream that only
    # greys the paper white (display Y 0.73 instead of ~0.95), so it is off.
    lightness_compression: bool = False
    # Remove the print's neutral-scale crossover: per-channel 1D remap of the
    # print dye densities (a crossover-free print stock), solved so every
    # grey of the scene prints neutral at unchanged luminance.
    neutralize: bool = False
    # Display-side colour intent, in Oklab: rotate the film's hue towards the
    # plain DRT's hue by this fraction (0 = pure film, 1 = DRT hues with the
    # film's lightness/chroma), and scale chroma.
    hue_preserve: float = 0.0
    chroma_gain: float = 1.0
    display: str = "sRGB - Display"
    view: str = "ACES 2.0 - SDR 100 nits (Rec.709)"


def _srgb_decode(enc: np.ndarray) -> np.ndarray:
    import colour
    return colour.cctf_decoding(np.fmax(enc, 0.0), function="sRGB")


def _srgb_to_oklab(enc: np.ndarray) -> np.ndarray:
    import colour
    return colour.XYZ_to_Oklab(colour.sRGB_to_XYZ(enc))


def _oklab_to_srgb(lab: np.ndarray) -> np.ndarray:
    import colour
    return colour.XYZ_to_sRGB(colour.Oklab_to_XYZ(lab))


class PrintDRTModel:
    def __init__(self, spec: PrintDRTSpec):
        import PyOpenColorIO as ocio
        from spektrafilm.profiles.io import load_profile
        from spektrafilm.runtime.params_builder import digest_params, init_params
        from spektrafilm.runtime.pipeline import SimulationPipeline

        self.spec = spec
        info = load_profile(spec.film_profile).info
        # Reversal (positive) film is viewed/scanned directly: no print stage.
        self.reversal = info.type == "positive"
        pr = spec.print_profile or info.target_print or "kodak_2383"
        p = init_params(film_profile=spec.film_profile, print_profile=pr)
        p.debug.lut_mode = True
        p.io.scan_film = self.reversal
        p.io.input_color_space = "ACES2065-1"
        p.io.input_cctf_decoding = False
        p.io.output_color_space = "sRGB"
        p.io.output_cctf_encoding = True
        if not spec.lightness_compression:
            from dataclasses import replace
            p.io.output_gamut_compress = replace(p.io.output_gamut_compress, lightness_compression=None)
        p = digest_params(p)
        # lut_mode switches the scanner references off; they only depend on
        # the stock reference densities (not on image content), so they are
        # safe inside a static LUT.
        p.scanner.white_correction = spec.white_correction
        p.scanner.black_correction = spec.black_correction
        p.scanner.white_level = spec.white_level
        p.scanner.black_level = spec.black_level
        self._pipeline = SimulationPipeline(p)
        cfg = ocio.Config.CreateFromBuiltinConfig(_OCIO_CONFIG)
        fwd = ocio.DisplayViewTransform(src="ACES2065-1", display=spec.display, view=spec.view)
        inv = ocio.DisplayViewTransform(src="ACES2065-1", display=spec.display, view=spec.view,
                                        direction=ocio.TRANSFORM_DIR_INVERSE)
        self._drt = cfg.getProcessor(fwd).getDefaultCPUProcessor()
        self._inv_drt = cfg.getProcessor(inv).getDefaultCPUProcessor()
        if self.reversal and spec.neutralize:
            raise ValueError("neutralize works on print densities; not available for reversal film")
        self._neutral_maps = self._solve_neutral_maps() if spec.neutralize else None

    # -- print density stage ------------------------------------------------

    def _cmy_print(self, aces: np.ndarray) -> np.ndarray:
        img = np.fmax(np.asarray(aces, dtype=np.float64).reshape(1, -1, 3), 0.0) * 2.0 ** self.spec.exposure_ev
        return np.asarray(self._pipeline.process(img, collect="cmy_print"), dtype=np.float64).reshape(-1, 3)

    def _scan(self, cmy: np.ndarray) -> np.ndarray:
        out = self._pipeline.process(np.asarray(cmy, dtype=np.float64).reshape(1, -1, 3),
                                     inject="cmy_print", collect="rgb_out")
        return np.asarray(out, dtype=np.float64).reshape(-1, 3)

    def _solve_neutral_maps(self, n: int = 241, iters: int = 15, h: float = 0.01):
        """Per-channel density maps D_c -> D'_c that print the grey scale neutral.

        For every grey of the scene, Newton-solve the three print dye densities
        whose scan is neutral (R = G = B) at the original luminance.
        """
        stops = np.linspace(-12.0, 12.0, n)
        grey = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
        d0 = self._cmy_print(grey)
        y0 = _srgb_decode(self._scan(d0)) @ np.array([0.2126, 0.7152, 0.0722])
        target = np.log(np.fmax(y0, 1e-6))[:, None]

        def resid(d):
            return np.log(np.fmax(_srgb_decode(self._scan(d)), 1e-6)) - target

        d = d0.copy()
        for _ in range(iters):
            f = resid(d)
            if np.max(np.abs(f)) < 1e-4:
                break
            jac = np.empty((n, 3, 3))
            for k in range(3):
                dk = d.copy()
                dk[:, k] += h
                jac[:, :, k] = (resid(dk) - f) / h
            step = np.linalg.solve(jac + 1e-6 * np.eye(3), f[..., None])[..., 0]
            d = d - np.clip(step, -0.2, 0.2)
        self.neutral_residual = np.abs(resid(d)).max(axis=1)
        maps = []
        for c in range(3):
            x, idx = np.unique(d0[:, c], return_index=True)
            maps.append((x, d[idx, c]))
        return maps

    def _apply_neutral_maps(self, cmy: np.ndarray) -> np.ndarray:
        out = np.empty_like(cmy)
        for c, (x, y) in enumerate(self._neutral_maps):
            v = cmy[:, c]
            r = np.interp(v, x, y)
            # outside the solved range keep the end offset (no flattening)
            r = np.where(v < x[0], v + (y[0] - x[0]), r)
            r = np.where(v > x[-1], v + (y[-1] - x[-1]), r)
            out[:, c] = r
        return out

    # -- display and LMT ----------------------------------------------------

    def display(self, aces: np.ndarray) -> np.ndarray:
        img = np.asarray(aces, dtype=np.float64)
        if self.reversal:
            rgb = self._pipeline.process(np.fmax(img.reshape(1, -1, 3), 0.0) * 2.0 ** self.spec.exposure_ev)
            out = np.clip(np.asarray(rgb, dtype=np.float64).reshape(-1, 3), 0.0, 1.0)
        else:
            cmy = self._cmy_print(img)
            if self._neutral_maps is not None:
                cmy = self._apply_neutral_maps(cmy)
            out = np.clip(self._scan(cmy), 0.0, 1.0)
        if self.spec.hue_preserve or self.spec.chroma_gain != 1.0:
            out = self._colour_intent(img.reshape(-1, 3), out)
        return out.reshape(img.shape)

    def _colour_intent(self, aces: np.ndarray, film: np.ndarray) -> np.ndarray:
        ref = np.clip(_apply_cpu(self._drt, aces * 2.0 ** self.spec.exposure_ev), 0.0, 1.0)
        lab_f, lab_r = _srgb_to_oklab(film), _srgb_to_oklab(ref)
        c_f = np.hypot(lab_f[:, 1], lab_f[:, 2])
        h_f = np.arctan2(lab_f[:, 2], lab_f[:, 1])
        c_r = np.hypot(lab_r[:, 1], lab_r[:, 2])
        h_r = np.arctan2(lab_r[:, 2], lab_r[:, 1])
        dh = (h_r - h_f + np.pi) % (2 * np.pi) - np.pi
        # hue of near-neutrals is undefined: fade the rotation in with chroma
        w = self.spec.hue_preserve * np.clip(np.minimum(c_r, c_f) / 0.03, 0.0, 1.0)
        h = h_f + w * dh
        c = c_f * self.spec.chroma_gain
        lab = np.stack([lab_f[:, 0], c * np.cos(h), c * np.sin(h)], axis=-1)
        return np.clip(_oklab_to_srgb(lab), 0.0, 1.0)

    def aces_out(self, aces: np.ndarray) -> np.ndarray:
        return _apply_cpu(self._inv_drt, self.display(aces))

    def acescc_lmt(self, acescc: np.ndarray) -> np.ndarray:
        return np.clip(aces_to_acescc(self.aces_out(acescc_to_aces(acescc))), ACESCC_MIN, ACESCC_MAX)
