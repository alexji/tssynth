from tssynth import marcs

model_atmosphere_file_1 = "/Users/alexji/lib/tssynth/tests/s5000_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
model_atmosphere_file_2 = "/Users/alexji/lib/tssynth/tests/p4000_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"
model_atmosphere_file_3 = "/Users/alexji/lib/tssynth/tests/p7750_g+4.5_m0.0_t02_st_z-0.25_a+0.10_c+0.00_n+0.00_o+0.10_r+0.00_s+0.00.mod"
model_atmosphere_file_4 = "/Users/alexji/lib/tssynth/tests/s7500_g+3.5_m1.0_t05_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"

def test_parse_marcs_model_1():
    header, model_structure = marcs.parse_marcs_model(model_atmosphere_file_1)
    assert header["Teff"] == 5000
    assert header["logg"] == 2.0
    assert header["vturb"] == 2.0
    assert header["feh"] == -2.0
    assert header["alphafe"] == 0.4
    assert header["spherical"] == True

    header, model_structure = marcs.parse_marcs_model(model_atmosphere_file_2)
    assert header["Teff"] == 4000
    assert header["logg"] == 4.5
    assert header["vturb"] == 1.0
    assert header["feh"] == 0.0
    assert header["alphafe"] == 0.0
    assert header["spherical"] == False

    header, model_structure, abundances, partial_pressures_lines = marcs.parse_marcs_model(model_atmosphere_file_1, get_all=True)
def test_parse_marcs_model_2():
    ## These ones have no spaces between adjacent columns because it's a formatted file!
    header, model_structure = marcs.parse_marcs_model(model_atmosphere_file_3)
    assert header["Teff"] == 7750
    assert header["logg"] == 4.5
    assert header["vturb"] == 2.0
    assert header["feh"] == -0.25
    assert header["alphafe"] == 0.1

    header, model_structure = marcs.parse_marcs_model(model_atmosphere_file_4)
    assert header["Teff"] == 7500
    assert header["logg"] == 3.5
    assert header["vturb"] == 5.0
    assert header["feh"] == -5.0
    assert header["alphafe"] == 0.4

def test_parse_marcs_model_3():
    """
    These are some edge cases that arose
    Have ****** or wrong numbers of columns
    """
    fnames = ["p3300_g+4.5_m0.0_t01_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
    "s2500_g+0.5_m1.0_t02_st_z-2.50_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
    "s6750_g+3.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
    "s7000_g+3.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
    "s2500_g+0.5_m1.0_t05_st_z-2.50_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
    "s5750_g+1.5_m1.0_t05_st_z-0.75_a+0.30_c+0.00_n+0.00_o+0.30_r+0.00_s+0.00.mod",
    "s8000_g+2.0_m1.0_t02_st_z-3.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"]
    for fname in fnames:
        header, model_structure = marcs.parse_marcs_model("/Users/alexji/lib/tssynth/tests/" + fname)
        

def test_interpolation():
    ## interpolate away from any grid points
    # Teff, logg, MH
    new_model_atmosphere_file = marcs.interpolate_marcs_model(5050, 2.1, -2.1)

    ## interpolate at a grid point
    new_model_atmosphere_file = marcs.interpolate_marcs_model(5000, 2.0, -2.0)
    