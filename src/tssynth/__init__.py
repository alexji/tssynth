import os, shutil

## HARDCODE during development
os.environ["TSEXEC_PATH"] = "/Users/alexji/lib/Turbospectrum_NLTE/exec-gf"
os.environ["TWD_BASE"] = "/Users/alexji/.tssynth"
os.environ["TSLINELIST_PATH"] = "/Users/alexji/lib/tssynth/data/linelists"
os.environ["TSDEPCOEFF_PATH"] = "/Users/alexji/bergemann_departure_coefficients"
os.environ["ALLMARCS_PATH"] = "/Users/alexji/MARCS/MARCS"

TSEXEC_PATH = os.environ.get('TSEXEC_PATH', None)
if TSEXEC_PATH is None:
    raise ValueError("Environment variable TSEXEC_PATH is not set.")
if not os.path.exists(TSEXEC_PATH):
    raise ValueError(f"{TSEXEC_PATH} does not exist.")
print(f"Using Turbospectrum from {TSEXEC_PATH}")

TWD_BASE = os.environ.get('TWD_BASE', None)
if TWD_BASE is None:
    raise ValueError("Environment variable TWD_BASE is not set.")
if not os.path.exists(TWD_BASE):
    print(f"Creating {TWD_BASE}")
    os.makedirs(TWD_BASE)
print(f"Using TWD_BASE at {TWD_BASE}")

TSLINELIST_PATH = os.environ.get('TSLINELIST_PATH', None)
if TSLINELIST_PATH is None:
    raise ValueError("Environment variable TSLINELIST_PATH is not set.")
if not os.path.exists(TSLINELIST_PATH):
    raise ValueError(f"{TSLINELIST_PATH} does not exist.")
print(f"Using linelists from {TSLINELIST_PATH}")

TSDEPCOEFF_PATH = os.environ.get('TSDEPCOEFF_PATH', None)
if TSDEPCOEFF_PATH is None:
    raise ValueError("Environment variable TSDEPCOEFF_PATH is not set.")
if not os.path.exists(TSDEPCOEFF_PATH):
    raise ValueError(f"{TSDEPCOEFF_PATH} does not exist.")
print(f"Using departure coefficients from {TSDEPCOEFF_PATH}")

ALLMARCS_PATH = os.environ.get('ALLMARCS_PATH', None)
if ALLMARCS_PATH is None:
    raise ValueError("Environment variable ALLMARCS_PATH is not set.")
if not os.path.exists(ALLMARCS_PATH):
    raise ValueError(f"{ALLMARCS_PATH} does not exist.")
print(f"Using MARCS atmospheres from {ALLMARCS_PATH}")

for executable in ['babsma_lu', 'bsyn_lu']:
    if not shutil.which(executable):
        raise EnvironmentError(f"{executable} not found in PATH. Please ensure Turbospectrum is installed and {executable} is in your PATH.")
print("Turbospectrum executables found in PATH.")

## Enable tssynth.solar_abundances and tssynth.solar_isotopes from tssynth
from .solar_abundances import solar_abundances, periodic_table
from .solar_isotopes import solar_isotopes

## Import basic interface
from .synthesizer import run_synth_lte
from . import utils
