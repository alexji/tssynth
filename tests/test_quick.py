import tssynth
from tssynth import synthesizer
import shutil

model_atmosphere_file = "/Users/alexji/lib/tssynth/tests/s5000_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"

def test_simple_methods():
    print(synthesizer.get_default_linelist_filenames())
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

def test_run_synth_lte_quick():
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

def test_run_synth_nlte_quick():
    # wmin, wmax, dw = 5090, 5100, 0.1
    # wave, norm, flux = synthesizer.run_synth_nlte(wmin, wmax, dw,
    #                           model_atmosphere_file=model_atmosphere_file)
    # assert len(wave) == len(norm) == len(flux) == round((wmax-wmin)/dw) + 1
    pass