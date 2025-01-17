from tssynth import marcs

model_atmosphere_file_1 = "/Users/alexji/lib/tssynth/tests/s5000_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
model_atmosphere_file_2 = "/Users/alexji/lib/tssynth/tests/p4000_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"

def test_header():
    header, abundances, model_structure, partial_pressures_lines = marcs.parse_marcs_model(model_atmosphere_file_1)
    assert header["Teff"] == 5000
    assert header["logg"] == 2.0
    assert header["vturb"] == 2.0
    assert header["feh"] == -2.0
    assert header["alphafe"] == 0.4
    assert header["spherical"] == True

    header, abundances, model_structure, partial_pressures_lines = marcs.parse_marcs_model(model_atmosphere_file_2)
    assert header["Teff"] == 4000
    assert header["logg"] == 4.5
    assert header["vturb"] == 1.0
    assert header["feh"] == 0.0
    assert header["alphafe"] == 0.0
    assert header["spherical"] == False
