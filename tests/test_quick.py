import tssynth
from tssynth import synthesizer
import os, shutil
import numpy as np

# import pytest
# @pytest.mark.quick
# @pytest.mark.full
# >> pytest -m fast

model_atmosphere_file = "/Users/alexji/lib/tssynth/tests/model_atmospheres/s5000_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"

def test_simple_methods():
    Teff, logg, vt, MH, aFe, spherical = synthesizer.parse_model_atmosphere_file_params(model_atmosphere_file)
    assert Teff == 5000
    assert logg == 2.0
    assert vt == 2.0
    assert MH == -2.0
    assert aFe == 0.4
    assert spherical == True
    linelistnames = synthesizer.get_default_linelist_filenames()
    assert len(linelistnames) == 14
    linelistnames = synthesizer.get_default_linelist_filenames(include_H=False)
    assert len(linelistnames) == 13

def test_run_synth_lte_quick_1():
    wmin, wmax, dw = 5090, 5100, 0.05
    twd = tssynth.utils.mkdtemp()
    wave, norm, flux = synthesizer.run_synth_lte(wmin, wmax, dw,
                              model_atmosphere_file=model_atmosphere_file,
                              twd=twd)
    assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
    
    XFedict = {"Fe": -0.2, "Mg": 0.5, "Ba": -1.2}
    wave, norm, flux = synthesizer.run_synth_lte(wmin, wmax, dw,
                              model_atmosphere_file=model_atmosphere_file,
                              XFedict=XFedict, twd=twd)
    assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
    shutil.rmtree(twd)

def test_run_synth_lte_quick_2():
    wmin, wmax, dw = 5090, 5100, 0.05
    Teff, logg, MH = 5050, 2.05, -2.05
    twd = tssynth.utils.mkdtemp()
    wave, norm, flux = synthesizer.run_synth_lte(wmin, wmax, dw,
                              Teff=Teff, logg=logg, MH=MH,
                              twd=twd)
    assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
    
    XFedict = {"Fe": -0.2, "Mg": 0.5, "Ba": -1.2}
    wave, norm, flux = synthesizer.run_synth_lte(wmin, wmax, dw,
                              Teff=Teff, logg=logg, MH=MH,
                              XFedict=XFedict, twd=twd)
    assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
    shutil.rmtree(twd)

def test_run_synth_lte_full():
    wmin, wmax, dw = 4000, 4500, 0.1
    twd = tssynth.utils.mkdtemp()
    wave, norm, flux = synthesizer.run_synth_lte(wmin, wmax, dw,
                                                Teff=4500, logg=1.0, MH=-1.5,
                                                twd=twd)
    fname = os.path.join(os.path.dirname(__file__), "spectra", "spectrum_0_LTE.spec")
    wavecomp, normcomp, fluxcomp = np.loadtxt(fname).T
    ii = (wavecomp >= wmin) & (wavecomp <= wmax)
    wavecomp, normcomp, fluxcomp = wavecomp[ii], normcomp[ii], fluxcomp[ii]
    check_wave = np.allclose(wave, wavecomp, atol=0.01)
    check_norm = np.allclose(norm, normcomp, atol=0.005)
    check_flux = np.allclose(flux, fluxcomp, rtol=0.005)
    errors = ""
    if not check_wave:
        errors += f"waves do not match\n"
    if not check_norm:
        diff = np.abs(norm-normcomp)
        errors += f"norms do not match, max err {diff.max()}, {np.sum(diff>0.005)}/{len(diff)} pixels failed\n"
    if not check_flux:
        ratio = np.abs(flux/fluxcomp)-1
        errors += f"fluxes do not match, max relerr {(np.abs(ratio).max())}, {np.sum(np.abs(ratio)>0.005)}/{len(ratio)} pixels failed\n"
    assert check_wave and check_norm and check_flux, errors
    shutil.rmtree(twd)

# def test_run_synth_nlte_quick():
#     wmin, wmax, dw = 5090, 5100, 0.1
#     wave, norm, flux = synthesizer.run_synth_nlte(wmin, wmax, dw,
#                               model_atmosphere_file=model_atmosphere_file)
#     assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
#     pass