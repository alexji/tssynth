# tssynth
Python wrapper for running turbospectrum (in LTE and NLTE)

Authors and Contributors
-------
 - Alex Ji (University of Chicago, main author)
 - TSFitPy Contributors (from which much of the code is based): Nick Storm, Jeffrey Gerber. Please see https://github.com/TSFitPy-developers/TSFitPy/ for citations.

## Installation
- Install https://github.com/bertrandplez/Turbospectrum_NLTE
  - Make sure you have `gfortran` in your path
  - Clone, go into `exec-gf`, and run `make`
  - There should now be `babsma_lu` and `bsyn_lu` in that directory
- Use tssynth_conda.yml to create a new conda environment and activate it
- Install tssynth
  - For now/development: clone the repository, `python -m pip install -e .`
  - Eventually: `pip install tssynth`
- Specify environment variables to needed paths (we specify this all in `src/tssynth/__init__.py` so they print out during imports; there may be a more elegant solution with objects to be done in the future)
  - `TSEXEC_PATH`: path to directory including `babsma_lu` and `bsyn_lu`
  - `TWD_BASE`: path to where temporary directories will be written
  - `TSLINELIST_PATH`: path to directory containing linelists (should be in `tssynth`)
  - `TSDEPCOEFF_PATH`: path to directory containing departure coefficients. Only needed for NLTE. This will be very big if you actually download everything.
  - `ALLMARCS_PATH`: path to directory containing full library of MARCS model atmospheres. Only needed if you will interpolate model atmospheres.
- Download relevant files (can use `tssynth.downloader`):
  - Linelists (Default VALD is included with `tssynth`)
  - MARCS model atmosphere grid
  - Departure coefficient grids for elements you want NLTE
- Compile the fortran interpolators
  - go to `tssynth/fortran` and run `make`
  - there should be `interpol_modeles` and `interpol_modeles_nlte`


## Usage (currently aspirational)
Method 1: specify stellar parameters

Use the MPIA Bergemann Group's default atomic data linelist to do this
```
import tssynth

Teff, logg, vt, MH = 5000, 2.00, 1.50, -2.00
wmin, wmax, dw = 5000, 5100, 0.1

# can specify XFedict as either element or atomic number
# Here we say [Mg/Fe] = +0.8, [Ti/Fe] = -0.5
XFedict = {"Mg":0.8, 22:-0.5}

# Synthesize LTE
wave, flux_lte = tssynth.run_synth_lte(wmin, wmax, dw, 
                                        Teff=Teff, logg=logg, vt=vt, MH=MH, 
                                        XFedict=XFedict)

# Synthesize NLTE
# This will check TSDEPCOEF_PATH for elements in NLTE_elements
wave, flux_nlte = tssynth.run_synth_nlte(wmin, wmax, dw, 
                                        Teff=Teff, logg=logg, vt=vt, MH=MH,
                                        XFedict=XFedict,
                                        NLTE_elements=["Mg","Fe"])
```

Method 2: specify stellar atmosphere file
```
## Interpolation will be figured out later
# atmosphere_fname = "atmosphere.txt"
# tssynth.interpolate_model_atmosphere(atmosphere_fname, Teff, logg, vt, MH)

wave, flux_lte = tssynth.run_synth_lte(wmin, wmax, dw,
                                        model_atmosphere="atmosphere.txt", 
                                        XFedict=XFedict)

wave, flux_nlte = tssynth.run_synth_nlte(wmin, wmax, dw,
                                        model_atmosphere="atmosphere.txt",
                                        XFedict=XFedict,
                                        NLTE_data_dir=NLTE_data_dir)
```

