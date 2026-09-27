# ACEScc → ACEScc negative-film LMT (ADX workflow)

The regular `spektrafilm-lut` bundles bake the full film → print → scan chain
into a **display-referred** cube. Applied as an LMT inside an ACES pipeline this
double-renders: the print + scan tone scale clips and flattens the image *before*
the DRT, so exposure and highlights are wrong.

`spektrafilm_lut_creator.aces_lmt` instead models a motion-picture film scan:

```
scene (ACES2065-1, "the real world")
  → camera exposure onto the virtual negative        spektrafilm filming stage (spectral upsampling, couplers)
  → development, CMY dye amounts                     spektrafilm `cmy_film` tap (no print stage)
  → spectral density of the processed negative       channel dye spectra + base/orange mask
  → Academy Printing Density                         SMPTE ST 2065-2 eq. (1), official responsivities
  → ADX16                                            SMPTE ST 2065-3: gains (1.00, 0.92, 0.95), Dmin → 1520
  → Academy ADX16 IDT                                OCIO builtin `ADX16_to_ACES2065-1`
  → ACES2065-1 → ACEScc (AP1)
```

The output is scene-referred: the negative's toe and shoulder are the only tone
compression, and the DRT (ACES 2.0 / OpenDRT) still does the display rendering.

## Calibration choices (all physical lab operations)

| step | default | what it models |
|---|---|---|
| `density_metric` | `apd` | `apd`: strict ST 2065-2. `print`: RP 180 printing density of spektrafilm's own printer (target print stock × neutral-filtered lamp), the metric the stock profiles were fitted against. Use `print` for still-photo negatives. |
| Dmin | always | scanner calibrated so the unexposed negative sits on the ADX aim (ADX16 1520 / ADX10 95). |
| `balance` | `grey` | Camera CC filtration: per-layer log-exposure offsets (Newton-solved) so an exposed 18 % grey card has equal CDD, then one scalar ACES exposure gain → grey = 0.18. Keeps Dmin on the aim, so the very steep IDT toe near Dmin is not disturbed. Alternatives: `printer` (printer-light density offsets; tints deep shadows), `aces` (per-channel gain after the IDT; drives saturated colours negative on badly balanced stocks), `none`. |
| `calibration` | `none` | `gamma`: per-stock scanner calibration (ST 2065-2 note 3), per-channel CDD gains equalising mid-scale gammas. Recommended for still-photo negatives, which the Academy IDT was not designed for. |

## Findings

* APD of Vision3 250D near grey has a slope of ≈0.53 per log E, and the Academy IDT assumes 0.55, so ADX and Vision3 agree closely.
  The printing density of 2383 × S_APD reproduces the official APD responsivities to within ~0.01 density, which validates the spectral integration.
* The profile dye spectra are NaN at 380–400 nm. Treating NaN as zero density made that band transparent and capped every
  printing density (a fake blue shoulder). The spectra are now edge-held (`_fill_edges`).
* Neutral-scale crossover (max−min channel, stops) of the baked 65³ cubes with alcedo-style trilinear sampling:

  | stock | −8 | −4 | −2 | +2 | +4 | +6 | +8 | +10 | peak out |
  |---|---|---|---|---|---|---|---|---|---|
  | Vision3 50D | .07 | .34 | .18 | .04 | .11 | .22 | .30 | .33 | 20 |
  | Vision3 200T | .13 | .26 | .14 | .04 | .09 | .19 | .24 | .26 | 17 |
  | Vision3 250D | .10 | .17 | .11 | .05 | .11 | .24 | .39 | .48 | 24 |
  | Vision3 500T | .10 | .14 | .09 | .06 | .14 | .29 | .41 | .49 | 24 |
  | Portra 400 (`print`,`gamma`) | .01 | .14 | .04 | .02 | .03 | .09 | .27 | .29 | 40 |

  Grey is preserved within 0.02 stop. The 65³ trilinear cube versus the exact model after the ACES 2.0 DRT on real raws:
  mean ΔE2000 0.18, p99 0.68.
* alcedo_studio ignores `DOMAIN_MIN/MAX` and clamps the ACEScc input to [0, 1]
  (linear 0.0012–223, i.e. −7.2…+10.3 stops around grey). Cubes are therefore baked on [0, 1].
  Below −7.2 stops the film is already on its base (toe floor ≈ −7.5…−8 stops), so little is lost.

## Scripts

```bash
python scripts/aces_lmt/bake.py kodak_vision3_250d kodak_vision3_500t            # → experiment_results/luts
python scripts/aces_lmt/bake.py kodak_portra_400 --metric print --calibration gamma
python scripts/aces_lmt/verify.py <cube> --film kodak_vision3_250d --raws <raws…>  # tone scale plot + ΔE
python scripts/aces_lmt/evaluate.py --auto-exposure --luts <cubes…> --raws <raws…>  # rawpy + OCIO ACES 2.0 contact sheets
python scripts/aces_lmt/probe_grey.py kodak_vision3_250d grey apd                    # stop-by-stop ADX table
```

The SMPTE ST 2065-2 data supplement (APD responsivities, influx spectrum) is downloaded from
pub.smpte.org on first use and cached in `~/.cache/spektrafilm/smpte_st2065_2` (override with
`SPEKTRAFILM_APD_DIR`). It is not redistributed in this repository.
