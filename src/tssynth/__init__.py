import os, shutil

## HARDCODE during development
os.environ["TURBOSPECTRUM_PATH"] = "/Users/alexji/lib/Turbospectrum_NLTE/exec-gf"
os.environ["TWD_BASE"] = "/Users/alexji/.tssynth"
os.environ["TSLINELIST_PATH"] = "/Users/alexji/lib/tssynth/data/linelists"

TURBOSPECTRUM_PATH = os.environ.get('TURBOSPECTRUM_PATH', None)
if TURBOSPECTRUM_PATH is None:
    raise ValueError("Environment variable TURBOSPECTRUM_PATH is not set.")
if not os.path.exists(TURBOSPECTRUM_PATH):
    raise ValueError(f"{TURBOSPECTRUM_PATH} does not exist.")
print(f"Using Turbospectrum from {TURBOSPECTRUM_PATH}")

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
