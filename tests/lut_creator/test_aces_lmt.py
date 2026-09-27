"""Tests for the ACEScc -> ACEScc negative-film (ADX) LMT."""
from __future__ import annotations

import numpy as np
import pytest

from spektrafilm_lut_creator import aces_lmt as lmt

pytest.importorskip("PyOpenColorIO")


def _apd_available() -> bool:
    try:
        lmt.load_apd_responsivities()
        return True
    except OSError:
        return False


needs_apd = pytest.mark.skipif(not _apd_available(), reason="SMPTE ST 2065-2 data not cached/reachable")


@pytest.mark.unit
def test_acescc_codec_matches_ocio_for_positive_values():
    rng = np.random.default_rng(0)
    aces = rng.lognormal(-2.0, 2.5, (4096, 3))
    aces = aces[(aces @ lmt._AP0_TO_AP1.T).min(axis=1) > 0]
    ocio = lmt._apply_cpu(lmt._ocio_cpu("ACEScc_to_ACES2065-1", inverse=True), aces)
    np.testing.assert_allclose(lmt.aces_to_acescc(aces), ocio, atol=2e-5)
    np.testing.assert_allclose(lmt.acescc_to_aces(lmt.aces_to_acescc(aces)), aces, rtol=1e-9)


@pytest.mark.unit
def test_acescc_negative_extension_round_trips():
    v = np.array([-0.5, -0.01, 0.0, 1e-6, 0.18, 100.0])
    np.testing.assert_allclose(lmt.acescc_decode(lmt.acescc_encode(v)), v, atol=1e-12)


@pytest.mark.unit
def test_adx16_idt_reference_point():
    # ACEScsc.Academy.ADX16_to_ACES REF_PT: neutral CDD 0.70 above Dmin -> 0.18.
    cv = np.full((1, 3), lmt.LAD_CDD * lmt.ADX16_SCALE + lmt.ADX16_OFFSET)
    np.testing.assert_allclose(lmt.adx16_to_aces(cv), 0.18, rtol=2e-3)


@pytest.mark.unit
def test_fill_edges_holds_nearest_valid_value():
    x = np.array([[np.nan, 1.0], [np.nan, 2.0], [3.0, np.nan], [4.0, 5.0]])
    np.testing.assert_array_equal(lmt._fill_edges(x), [[3.0, 1.0], [3.0, 2.0], [3.0, 3.5], [4.0, 5.0]])


@needs_apd
@pytest.mark.unit
def test_printing_density_of_neutral_filter_is_its_density():
    wl, pi = lmt.load_apd_responsivities()
    film_wl = np.arange(380.0, 781.0, 5.0)
    dens = np.full((2, film_wl.size), 0.8)
    np.testing.assert_allclose(lmt.printing_density(dens, film_wl, wl, pi), 0.8, atol=1e-9)


@needs_apd
@pytest.mark.integration
@pytest.mark.parametrize("balance", ["grey", "printer", "aces"])
def test_negative_model_preserves_mid_grey(balance):
    model = lmt.NegativeADXModel(lmt.NegativeADXSpec(film_profile="kodak_vision3_250d", balance=balance))
    grey = np.full((1, 3), 0.18)
    np.testing.assert_allclose(model.aces_out(grey), grey, rtol=1e-3)
    # Dmin (unexposed negative) sits on the ADX16 aim.
    np.testing.assert_allclose(model.adx16(np.zeros((1, 3))), (lmt.ADX16_OFFSET + model.cdd_offset * lmt.ADX16_SCALE)[None], atol=1e-6)


@needs_apd
@pytest.mark.integration
def test_negative_model_keeps_scene_dynamic_range():
    model = lmt.NegativeADXModel(lmt.NegativeADXSpec(film_profile="kodak_vision3_250d"))
    stops = np.array([-4.0, -2.0, 2.0, 4.0, 6.0])
    aces = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
    out = np.log2((model.aces_out(aces) @ lmt._AP0_TO_AP1.T) / 0.18)
    # Near-unity tone scale around grey, highlights above display white survive.
    np.testing.assert_allclose(out.mean(axis=1)[:4], stops[:4], atol=0.35)
    assert out.mean(axis=1)[-1] > 5.0
    assert np.ptp(out, axis=1).max() < 0.35  # neutral scale stays neutral


@needs_apd
@pytest.mark.integration
def test_bake_shape_and_cube_writer(tmp_path):
    model = lmt.NegativeADXModel(lmt.NegativeADXSpec(film_profile="kodak_vision3_250d"))
    table = lmt.bake_acescc_lmt(model, size=5, domain=(0.0, 1.0))
    assert table.shape == (5, 5, 5, 3)
    # neutral diagonal stays neutral
    diag = table[np.arange(5), np.arange(5), np.arange(5)]
    assert np.ptp(diag[1:4], axis=1).max() < 0.02
    path = tmp_path / "t.cube"
    lmt.write_cube(table, path, title="t", domain=(0.0, 1.0))
    text = path.read_text()
    assert "LUT_3D_SIZE 5" in text and "DOMAIN_MAX 1 1 1" in text
    # SPEKTRAFILM_LICENSE.txt: every LUT names the author, source and license
    assert "Andrea Volpato" in text and "github.com/andreavolpato/spektrafilm" in text and "CC BY-SA 4.0" in text
    lmt.write_lut_license_files(tmp_path)
    assert (tmp_path / "SPEKTRAFILM_LICENSE.txt").exists() and (tmp_path / "CHANGELOG.txt").exists()


@pytest.mark.integration
@pytest.mark.parametrize("film", ["kodak_vision3_250d", "fujifilm_velvia_100"])
def test_print_drt_neutralize_makes_grey_scale_neutral(film):
    import colour
    model = lmt.PrintDRTModel(lmt.PrintDRTSpec(film_profile=film, neutralize=True))
    base = lmt.PrintDRTModel(lmt.PrintDRTSpec(film_profile=film))
    stops = np.array([-4.0, -2.0, 0.0, 2.0, 4.0])
    grey = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
    lab = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(model.display(grey)))
    lab0 = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(base.display(grey)))
    assert np.abs(lab[:, 1:]).max() < 0.3
    np.testing.assert_allclose(lab[:, 0], lab0[:, 0], atol=0.2)  # lightness unchanged
    # inverse DRT makes it scene-referred: forward DRT reproduces the display
    np.testing.assert_allclose(lmt._apply_cpu(model._drt, model.aces_out(grey)), model._to_display_code(model.display(grey)), atol=2e-3)
