# spektrafilm-lut — ACES LMT export for spektrafilm

> [!IMPORTANT]
> **This repository is a fork of [spektrafilm](https://github.com/andreavolpato/spektrafilm) by Andrea Volpato.**
>
> - spektrafilm is the work of **Andrea Volpato**. All physical models, film profiles, print profiles and the runtime come from his project.
> - The source code is licensed under the **GNU GPLv3** ([LICENSE](LICENSE)).
> - The film profiles and all LUTs made from them are licensed under **CC BY-SA 4.0** with the spektrafilm preamble ([SPEKTRAFILM_LICENSE.txt](SPEKTRAFILM_LICENSE.txt)). Each LUT must name Andrea Volpato and link to https://github.com/andreavolpato/spektrafilm.
> - To cite spektrafilm, use [CITATION.cff](CITATION.cff).
> - This fork is not an official spektrafilm release. Andrea Volpato does not maintain or endorse it. Send questions about this fork to this repository, not to the original project.
> - You can support the original author here: [Buy me a coffee](https://buymeacoffee.com/andreavolpato).
>
> **Original README:** [README_UPSTREAM.md](README_UPSTREAM.md). Read it for the spektrafilm model, the GUI and the full documentation.

This fork adds one function to spektrafilm: it exports the film simulation as an **ACEScc → ACEScc look modification transform (LMT)**. You apply this LMT as a 3D LUT in an ACES pipeline, before the output transform (DRT). The target application is [Alcedo Studio](https://github.com/zidage/AlcedoStudio).

## Where the changes are

| Branch | Contents |
|---|---|
| `aces-adx-lmt` | **All changes of this fork.** Use this branch. |
| `main` | The upstream spektrafilm code, plus this README and the README backup. There are no code changes. |

The fork starts at upstream commit `3bb2c2d` ("fix: vlog midgray exposure"). The `upstream` remote points to https://github.com/andreavolpato/spektrafilm.

## Background: the problem

Alcedo Studio applies a look LUT in ACEScc (AP1 primaries), before the DRT. The DRT is ACES 2.0 or OpenDRT.

The spektrafilm LUT creator makes **display-referred** LUTs. Such a LUT contains the full film → print → scan chain, so its output is already a finished display image. When you apply it as an LMT, the DRT then does a second tone mapping on that image. We measured this result:

- The paper white of the print stays at display L\* 84–88. The highlights never become white.
- The tone scale becomes flat above +3 stops. The image looks clipped and the exposure looks incorrect.

This fork gives two solutions. Section 1 is a technical film scan. Section 2 is the recommended look.

## 1. ADX film-scan LMT (technical scan)

Code: `NegativeADXModel` in [src/spektrafilm_lut_creator/aces_lmt.py](https://github.com/zidage/spektrafilm-lut/blob/aces-adx-lmt/src/spektrafilm_lut_creator/aces_lmt.py).

This model follows the motion-picture scan workflow. It uses the negative only. There is no print stage.

```
ACEScc → ACES2065-1 (the scene, the "real world")
  → spektrafilm: expose the negative, develop it (tap cmy_film)
  → spectral density of the negative (dye spectra + base/mask)
  → Academy Printing Density, APD (SMPTE ST 2065-2, equation 1)
  → ADX16 code values (SMPTE ST 2065-3: gains 1.00/0.92/0.95, Dmin at 1520)
  → Academy ADX16 IDT (OCIO builtin ADX16_to_ACES2065-1)
  → ACES2065-1 → ACEScc
```

Calibration steps. Each step is a real laboratory operation.

- **Dmin:** The virtual scanner puts the unexposed negative at the ADX aim (ADX16 1520).
- **Grey balance (`balance="grey"`):** Per-layer camera CC filtration makes an exposed 18 % grey card neutral. A Newton solver finds the filtration. Then one scalar ACES gain puts the grey at 0.18.
- **Printing-density metric:** `apd` uses the official ST 2065-2 responsivities. `print` uses the target print stock of spektrafilm under its printer lamp.
- **Per-stock calibration (`calibration="gamma"`):** This is optional. It equalizes the mid-scale channel gammas. Use it for still-photo negatives.

Results:

- The printing-density gamma of Vision3 250D is approximately 0.53. The Academy IDT expects 0.55. Thus ADX and Vision3 agree.
- Mid grey stays at 0.18 (error < 0.02 stop).
- The negative keeps the scene dynamic range: +10 stops of scene exposure give approximately +7 stops of output, with a smooth film shoulder.
- **Limit:** The IDT makes the negative tone scale linear again. The result is a flat technical image with a mid-scale contrast of 1.10, which is lower than ACES 2.0 alone (1.19). The "film look" comes from the print, and this model has no print.

The first run downloads the SMPTE ST 2065-2 data supplement from pub.smpte.org. The code caches it in `~/.cache/spektrafilm/smpte_st2065_2`. You can change this folder with the `SPEKTRAFILM_APD_DIR` environment variable. This repository does not contain the SMPTE data.

## 2. Print-chain LMT through the inverse ACES 2.0 output transform (recommended)

Code: `PrintDRTModel` in [src/spektrafilm_lut_creator/aces_lmt.py](https://github.com/zidage/spektrafilm-lut/blob/aces-adx-lmt/src/spektrafilm_lut_creator/aces_lmt.py).

```
LMT = (ACES 2.0 SDR output transform)⁻¹ ∘ scan( print( negative(scene) ) )
```

In Alcedo Studio, the ACES 2.0 output transform comes after the LMT. The inverse in the LUT and the forward transform in Alcedo Studio cancel each other. Thus the screen shows the spektrafilm print image. The LMT only encodes this image as scene-referred ACEScc data. Camera "film simulation" LUTs that are converted in DaVinci Resolve use the same structure.

Changes relative to the spektrafilm GUI rendering:

- **Scanner white and black references are on.** The paper white goes to display 0.98. The print Dmax goes to display 0.005. These references come from the film and print reference densities only, not from image content, so a static LUT can contain them.
- **The runtime lightness roll-off is off** (`lightness_compression` of the cam16ucs output gamut compression). With this roll-off on, the paper white is at display Y 0.73 (L\* 88) and the highlights look dull. With it off, the paper white is at Y 0.95 (L\* 98).
- **Reversal film** (Velvia, Provia, Ektachrome, Kodachrome) does not go through a print. The model scans the slide directly.
- **Display:** The inverse is made for ACES 2.0 SDR 100 nit, Rec.709 primaries, gamma 2.2 encoding. This is the default display setting of Alcedo Studio.

Measured results after ACES 2.0 SDR (script `analyze_luts.py`):

| LMT | Mid grey L\* | System contrast at grey | L\* at +4 / +6 stops |
|---|---|---|---|
| No LMT (ACES 2.0 only) | 37.8 | 1.19 | 88 / 97 |
| Fujifilm film simulations, converted in Resolve (Pro Neg, Provia, Eterna) | 45.7–45.9 | 1.44–1.45 | 91–92 / 99 |
| ADX film-scan LMT, Vision3 250D | 37.8 | 1.10 | 87 / 95 |
| Print-chain LMT, Vision3 250D → 2383 | 48.1 | 1.32 | 95 / 98 |
| Print-chain LMT, 36 film/print pairs | 46–49 | 1.25–2.07 | 95–98 / 97–98 |

The 65³ cube with trilinear sampling agrees with the exact model: mean ΔE2000 0.07, p99 0.46, in the valid input range.

### Optional print corrections (off by default)

These options exist in `PrintDRTSpec`. They are off by default because they remove part of the film character.

- **`neutralize`:** For reversal film, the solver works on the dye densities of the slide. The 2383 print has a crossover: the shadows at −2 stops are green (a\* −4.7) and the highlights at +2 to +4 stops are yellow (b\* +5 to +6.5). This option adds one 1D curve per channel on the print dye densities. A Newton solver makes each scene grey print neutral at the same L\*. The result is a crossover-free print stock.
- **`hue_preserve`** (0 to 1): This option turns the film hue in Oklab toward the hue of ACES 2.0 alone. It keeps the film lightness and chroma. At 0.5, the mean hue error on a ColorChecker goes from 8.9° to 4.2°.
- **`chroma_gain`:** This option multiplies the Oklab chroma.

## Files that this fork adds or changes

| File | Purpose |
|---|---|
| `src/spektrafilm_lut_creator/aces_lmt.py` | New. Both models, APD/ADX math, ACEScc codec, LUT sampling, `.cube` writer with license header. |
| `scripts/aces_lmt/bake_print_drt.py` | New. Makes print-chain LMTs (section 2). |
| `scripts/aces_lmt/bake_alcedo_set.sh` | New. Makes the full set of 36 print-chain LMTs for Alcedo Studio. |
| `scripts/aces_lmt/bake.py` | New. Makes ADX film-scan LMTs (section 1). |
| `scripts/aces_lmt/evaluate.py` | New. Renders camera raw files (rawpy) through the LMTs and ACES 2.0 (OCIO), with the LUT sampling of Alcedo Studio. |
| `scripts/aces_lmt/verify.py`, `analyze_luts.py`, `analyze_print_chain.py`, `probe_*.py` | New. Tone-scale, colour and ΔE measurements. |
| `scripts/aces_lmt/README.md` | New. Technical notes and measurements. |
| `tests/lut_creator/test_aces_lmt.py` | New. Unit and integration tests. |
| `.gitignore` | Changed. Ignores `/experiment_results/`. |
| `README.md`, `README_UPSTREAM.md` | This README. The upstream README is kept unchanged as `README_UPSTREAM.md`. |

The fork does not change the spektrafilm runtime, the film profiles or the print profiles.

## Get the code and set it up

You need Git, Python 3.13 and [uv](https://docs.astral.sh/uv/). The commands are for Git Bash on Windows. On Linux and macOS, use `.venv/bin/python` instead of `.venv/Scripts/python.exe`.

1. Clone the repository:

   ```bash
   git clone https://github.com/zidage/spektrafilm-lut.git
   ```

2. Go into the folder:

   ```bash
   cd spektrafilm-lut
   ```

3. Change to the branch with the changes:

   ```bash
   git checkout aces-adx-lmt
   ```

4. Make a Python 3.13 virtual environment:

   ```bash
   uv venv --python 3.13 .venv
   ```

5. Install spektrafilm with the development tools. This also installs OpenColorIO:

   ```bash
   uv pip install --python .venv/Scripts/python.exe -e ".[dev]"
   ```

6. Run the tests. The first run downloads the SMPTE data, so you need an internet connection:

   ```bash
   .venv/Scripts/python.exe -m pytest tests/lut_creator/test_aces_lmt.py -q
   ```

To get upstream updates, add the original repository as a remote and fetch it:

```bash
git remote add upstream https://github.com/andreavolpato/spektrafilm.git
```

```bash
git fetch upstream
```

## Make the LUTs

The scripts write the LUTs to `experiment_results/`. Git ignores this folder. The repository does not contain `.cube` files. You make them on your computer.

### File names

The scripts give each print-chain LMT a short name: `Brand_Film_Print[_calib].cube`, for example `Kodak_Vision3-250D_2383_NH.cube`. The LUT title in the file is the same name.

| Part | Values |
|---|---|
| Brand | `Kodak`, `Fuji` |
| Film | `Vision3-250D`, `Portra400`, `Portra800P1` (push 1), `Pro400H`, `Velvia100`, ... |
| Print | `2383`, `2393` (Kodak Vision print film), `Endura` (Kodak Portra Endura paper), `CA` (Fujifilm Crystal Archive Type II paper), `Slide` (reversal film, no print) |
| Calib | No tag: no correction. `N`: `neutralize`. `H`: `hue_preserve` (the value is in the file header). `C`: `chroma_gain`. |

Use `--long-names` to get the old descriptive file names.

### Make the full set for Alcedo Studio

This script makes 36 LUTs: 5 cine negatives × 2 print films, 8 Kodak and 3 Fujifilm still negatives × 2 papers, and 4 reversal films. Give the output folder and the options:

```bash
scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/plain
```

```bash
scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/N --neutralize
```

```bash
scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/NH --neutralize --hue-preserve 0.25
```

The `NH` set uses a small hue correction (0.25). It keeps most of the film colour. On a ColorChecker, the mean hue difference to ACES 2.0 changes from 8.9° to 6.1° for Vision3 250D → 2383.

### Make one LUT

Make one print-chain LMT (recommended):

```bash
.venv/Scripts/python.exe scripts/aces_lmt/bake_print_drt.py kodak_vision3_250d --out experiment_results/luts_final
```

Select a different print or paper with `--print`:

```bash
.venv/Scripts/python.exe scripts/aces_lmt/bake_print_drt.py kodak_portra_400 --print fujifilm_crystal_archive_typeii --out experiment_results/luts_final
```

Make a reversal-film LMT. The script finds the film type automatically:

```bash
.venv/Scripts/python.exe scripts/aces_lmt/bake_print_drt.py fujifilm_velvia_100 --out experiment_results/luts_final
```

Make an ADX film-scan LMT (technical scan):

```bash
.venv/Scripts/python.exe scripts/aces_lmt/bake.py kodak_vision3_250d --out experiment_results/luts
```

Each run also writes `SPEKTRAFILM_LICENSE.txt` and `CHANGELOG.txt` into the output folder. Each `.cube` file has this header:

```
# Derived from spektrafilm by Andrea Volpato
# https://github.com/andreavolpato/spektrafilm
# Licensed CC BY-SA 4.0 (see SPEKTRAFILM_LICENSE.txt)
# Modified by zidage: ACES LMT export (https://github.com/zidage/spektrafilm-lut, branch aces-adx-lmt)
```

Do not remove this header or these two files when you share LUTs. The spektrafilm LUT license requires them.

The film and print names are the file names in `src/spektrafilm/data/profiles/`, without `.json`. One LUT takes approximately 3 seconds.

## Use the LUTs in Alcedo Studio

[Alcedo Studio](https://github.com/zidage/AlcedoStudio) applies a LUT in a Color Grade node, in ACEScc. Do these steps:

1. Make the LUTs (see the section above).
2. In Alcedo Studio, open the LUT panel of a Color Grade node.
3. Click **Open LUT folder**.
4. Copy the `.cube` files into this folder. Also copy `SPEKTRAFILM_LICENSE.txt` and `CHANGELOG.txt`.
5. Click **Refresh LUT catalog**.
6. Select the LUT in the Color Grade node.
7. In **Display Transform**, select **ACES 2.0**. Keep the default display: Rec.709 primaries, **Gamma 2.2** encoding, peak luminance 100 nits.
8. Adjust the **Exposure** control so that a mid-grey subject is near 0.18 scene-linear. Alcedo Studio sets the sensor clip to 1.0 and does not set mid grey automatically. The LUT expects a correctly exposed scene, as a film camera does.

Rules for correct results:

- Use print-chain LMTs (`*_invACES2_*`) **only with ACES 2.0 SDR, Rec.709, gamma 2.2, 100 nits**. With OpenDRT, HDR or another encoding, the inverse does not match and the image is incorrect.
- ADX film-scan LMTs (`*_adx_*`) do not contain a DRT inverse. You can use them with each DRT. They give a technical scan without a print look.
- Alcedo Studio uses ACEScc [0, 1] as the LUT input range. This is scene-linear 0.0012 to 223 (−7.2 to +10.3 stops from mid grey). It clamps values outside this range. The scripts sample the LUTs in this range.
- Alcedo Studio uses trilinear sampling. `evaluate.py` uses the same sampling. Thus the test renders show the result of Alcedo Studio.
- In a Color Grade node, the controls before the LUT (white balance, exposure, contrast, curves, saturation) change the "scene" before the film. Shadows and Highlights come after the LUT.

## Test the LUTs on camera raw files

This command renders raw files through ACES 2.0 with and without the LMTs. It saves contact sheets in `experiment_results/renders/`:

```bash
.venv/Scripts/python.exe scripts/aces_lmt/evaluate.py --auto-exposure --luts experiment_results/luts_final/*.cube --raws path/to/raw1.CR3 path/to/raw2.NEF
```

This command measures the tone scale, the system contrast and the ColorChecker colours of LUTs:

```bash
.venv/Scripts/python.exe scripts/aces_lmt/analyze_luts.py --luts experiment_results/luts_final/*.cube
```

## Known limits

- **The print-chain LMTs are tied to one output transform** (ACES 2.0 SDR 100 nit, Rec.709, gamma 2.2). For another display, make new LUTs with a different `display` in `PrintDRTSpec`.
- **The print sets the dynamic range.** Highlights go smoothly to white at the print shoulder, near +6 stops. This agrees with a real print. For the full negative range, use the ADX film-scan LMT.
- **Some colour differences are model errors.** The profiles use generic dye spectra. This causes hue shifts, for example the blue sky turns approximately −10° to −15° toward cyan. Some of these shifts are real film character, and some are model errors. The options `neutralize` and `hue_preserve` can correct them.
- **Still-photo negatives in the ADX model** have 0.5–1 stop crossovers, because the Academy IDT is made for motion-picture negatives. Use `--metric print --calibration gamma` for these films.

## License of this fork

- New and changed source code: GNU GPLv3, the same as spektrafilm.
- LUTs that you make with this code: CC BY-SA 4.0 under [SPEKTRAFILM_LICENSE.txt](SPEKTRAFILM_LICENSE.txt). Keep the attribution header, the license file and `CHANGELOG.txt` with each copy. Do not sell the LUTs as a LUT pack. This is the wish of the original author.
- SMPTE ST 2065-2 data: the code downloads it from pub.smpte.org and does not redistribute it.

Original project: **spektrafilm by Andrea Volpato**, https://github.com/andreavolpato/spektrafilm.
