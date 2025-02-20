import numpy as np
import os, sys, shutil
import subprocess
from . import utils, marcs
from .solar_abundances import solar_abundances_Z

TSEXEC_PATH = os.environ.get('TSEXEC_PATH', None)
TSDATA_PATH = os.environ.get('TSDATA_PATH', os.path.join(TSEXEC_PATH, '../DATA'))
TWD_BASE = os.environ.get('TWD_BASE', None)
TSLINELIST_PATH = os.environ.get('TSLINELIST_PATH', None)
TSDEPCOEFF_PATH = os.environ.get('TSDEPCOEFF_PATH', None)
ALLMARCS_PATH = os.environ.get('ALLMARCS_PATH', None)

def run_synth_lte(wmin, wmax, dw,
                  Teff=None, logg=None, vt=2.0, MH=None, aFe=None,
                  model_atmosphere_file=None,
                  linelist_filenames=None,
                  XFedict=None, 
                  modelopac_file=None,
                  twd=None, delete_twd=False,
                  spherical=None, verbose=False):
    """
    Run LTE spectrum synthesis with Turbospectrum.

    Parameters:
    wmin (float): Minimum wavelength (A)
    wmax (float): Maximum wavelength (A)
    dw (float): Wavelength step (A)
    Teff (float): Effective temperature (default: None)
    logg (float): Surface gravity (default: None)
    vt (float): Microturbulent velocity (default: None)
    MH (float): Metallicity (default: None)
    aFe (float): Alpha enhancement (default: None)
    model_atmosphere_file (str): Path to the model atmosphere file (default: None)
    linelist_filenames (list or str): List of line list filenames or a single filename
        (default: None, uses get_default_linelist_filenames in TSLINELIST_PATH)
    XFedict (dict): Dictionary of element abundances (default: None)
        Turbospectrum assumes [X/Fe] following the model atmosphere if not specified.
    modelopac_file (str): Path to the model opacity file (default: None)
        This is useful if you want to compute many spectra from one model atmosphere/composition.
    twd (str): Temporary working directory (default: None, creates a new one with utils.mkdtemp)
        All the work is done in this directory.
    delete_twd (bool): Delete the temporary working directory after the function finishes (default: False)

    Turbospectrum runs in two steps.
    (1) babsma_lu: Computes the model opacity
    (2) bsyn_lu: Computes the spectrum
    It needs a temporary directory to work in.
    
    Returns:
    tuple: wave (numpy array), norm (numpy array), flux (numpy array)
    """

    if twd is None:
        twd = utils.mkdtemp()
        sys.stdout.write(f"Temporary working directory: {twd}")

    ## Model Atmosphere File
    if any(param is not None for param in [Teff, logg, MH]):
        assert vt == 2.0, f"for now vt must be 2.0 for atmosphere interpolation, specified {vt}"
        ## Specify the parameters here
        if not all(param is not None for param in [Teff, logg, MH]):
            raise ValueError("If any of Teff, logg, or MH is provided, they all need to be provided.")
        if aFe is None:
            if MH < -1.0: aFe = 0.4
            elif -1.0 < MH < 0.0: aFe = -0.4 * MH
            else: aFe = 0.0
        ## Automatically choose between spherical and plane parallel, arbitrary bound for now
        if spherical is None: spherical = logg <= 3.25
        elif spherical:
            assert logg <= 3.5, f"for spherical models, logg <= 3.5, specified {logg}"
        else:
            assert logg >= 3.5, f"for plane parallel models, logg >= 3.5, specified {logg}"
        
        ## Interpolate a model atmosphere and save into the twd
        model_atmosphere_file = os.path.join(twd, "marcs.interp")
        marcs.interpolate_marcs_model(Teff, logg, MH, model_atmosphere_file, spherical=spherical)
        is_marcsfile=False
    else:
        ## Specify a model atmosphere file
        assert model_atmosphere_file is not None, model_atmosphere_file
        assert os.path.exists(model_atmosphere_file), model_atmosphere_file
        Teff, logg, vt, MH, aFe, spherical = parse_model_atmosphere_file_params(model_atmosphere_file)
        is_marcsfile=True
    
    ## Line List
    if linelist_filenames is None:
        linelist_filenames = get_default_linelist_filenames()
    elif isinstance(linelist_filenames, str):
        linelist_filenames = [linelist_filenames]
    # check they're all here
    missing_files = [filename for filename in linelist_filenames if not os.path.exists(filename)]
    if missing_files:
        raise FileNotFoundError(f"Line list files not found: {', '.join(missing_files)}")
    
    ## Individual Abundances
    if XFedict is None:
        indiv_abu = {}
    else:
        indiv_abu = utils.parse_XFe_dict(XFedict)
    
    ## Set up the working directory
    if not os.path.exists(os.path.join(twd, 'DATA')):
        os.symlink(TSDATA_PATH, os.path.join(twd, 'DATA'))

    ## Run babsma_lu for Model Opacity
    kws_babsma_lu = dict(twd=twd,
                         wmin=wmin, wmax=wmax, dwl=dw,
                         modelfilename=model_atmosphere_file,
                         modelopacname=None,
                         MH=MH, aFe=aFe, indiv_abu=indiv_abu,
                         vt=vt, spherical=spherical,
                         is_marcsfile=is_marcsfile,
                         verbose=verbose)
    if modelopac_file is None:
        try:
            modelopac_file = run_babsma_lu(**kws_babsma_lu)
        except Exception as e:
            print("twd", twd)
            print(e)
            raise
    else: ## This is if you want to compute many spectra from one model opacity
        assert os.path.exists(modelopac_file), modelopac_file
        sys.stdout.write(f"Using model opacity file: {modelopac_file}, copying to twd")
        shutil.copy(modelopac_file,twd)
        modelopac_file= os.path.join(twd,os.path.basename(modelopac_file))
    
    ## Run bsyn_lu for Spectrum
    kws_bsyn_lu = kws_babsma_lu.copy()
    kws_bsyn_lu["modelopacname"] = modelopac_file
    kws_bsyn_lu["costheta"] = 1.0
    kws_bsyn_lu["isotopes"] = {}
    kws_bsyn_lu["linelistfilenames"] = linelist_filenames

    outfilename = run_bsyn_lu(**kws_bsyn_lu)

    ## Read the spectrum output
    turboOut= np.loadtxt(outfilename)
    if len(turboOut) == 0:
        raise RuntimeError("Turbospectrum output is empty: something did not work")

    wave = turboOut[:,0]
    norm = turboOut[:,1]
    flux = turboOut[:,2]


    return wave, norm, flux

def run_babsma_lu(twd, wmin, wmax, dwl,
                  modelfilename, modelopacname,
                  MH, aFe, indiv_abu, vt, spherical,
                  is_marcsfile=True,
                  verbose=False):
    """
    - create babsma_lu parameter file
    - call basbma_lu from TSEXEC_PATH
    - return filename of modelopac
    """
    scriptfilename= os.path.join(twd,'babsma.par')
    modelopacname= os.path.join(twd,'mopac')
    _write_script(scriptfilename,
                    wmin,wmax,dwl,
                    None,
                    modelfilename,
                    is_marcsfile,
                    modelopacname,
                    MH,
                    aFe,
                    indiv_abu,
                    vt,
                    spherical,
                    None,None,None,bsyn=False)
    # Run babsma
    sys.stdout.write('\r'+"Running Turbospectrum babsma_lu ...\r")
    sys.stdout.flush()
    if verbose:
        stdout= None
        stderr= None
    else:
        stdout= open('/dev/null', 'w')
        stderr= subprocess.STDOUT
    try:
        p= subprocess.Popen([os.path.join(TSEXEC_PATH, 'babsma_lu')],
                            cwd=twd,
                            stdin=subprocess.PIPE,
                            stdout=stdout,
                            stderr=stderr)
        with open(os.path.join(twd,'babsma.par'),'r') as parfile:
            for line in parfile:
                p.stdin.write(line.encode('utf-8'))
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        #for linelistfilename in linelistfilenames:
        #    os.remove(linelistfilename,twd)
        #if os.path.exists(os.path.join(twd,'DATA')):
        #    os.remove(os.path.join(twd,'DATA'))
        raise RuntimeError("Running babsma_lu failed ...")
    finally:
        #if os.path.exists(os.path.join(twd,'babsma.par')) \
        #   and outfname is None: #not 'saveTurboInput' in kwargs:
        #    os.remove(os.path.join(twd,'babsma.par'))
        # sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
    #if isinstance(modelopac,str):
    #    shutil.copy(modelopacname,modelopac)
    return modelopacname

def run_bsyn_lu(twd, wmin, wmax, dwl, costheta, modelfilename, is_marcsfile, 
                modelopacname, MH, aFe, indiv_abu, vt, spherical,
                isotopes, linelistfilenames, outfname=None, verbose=False):
    """
    """
    scriptfilename= os.path.join(twd,'bsyn.par')
    outfilename= os.path.join(twd,'bsyn.out')
    _write_script(scriptfilename,
                  wmin,wmax,dwl,
                  costheta,
                  modelfilename,
                  is_marcsfile,
                  modelopacname,
                  MH,
                  aFe,
                  indiv_abu,
                  vt,
                  spherical,
                  outfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=True)
    # Run bsyn
    sys.stdout.write('\r'+"Running Turbospectrum bsyn_lu ...\r")
    sys.stdout.flush()
    if verbose:
        stdout= None
        stderr= None
    else:
        stdout= open('/dev/null', 'w')
        stderr= subprocess.STDOUT
    try:
        p= subprocess.Popen([os.path.join(TSEXEC_PATH, 'bsyn_lu')],
                            cwd=twd,
                            stdin=subprocess.PIPE,
                            stdout=stdout,
                            stderr=stderr)
        with open(os.path.join(twd,'bsyn.par'),'r') as parfile:
            for line in parfile:
                p.stdin.write(line.encode('utf-8'))
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        raise RuntimeError("Running bsyn_lu failed ...")
    finally:
        if outfname is not None:
            turbosavefilename= outfname
            if os.path.dirname(turbosavefilename) == '':
                turbosavefilename= os.path.join(os.getcwd(),turbosavefilename)
            try:
                subprocess.check_call(['tar','cvzf',turbosavefilename,
                                       os.path.basename(os.path.normpath(twd))])
            except subprocess.CalledProcessError:
                raise RuntimeError("Tar-zipping the Turbospectrum input and output failed; you will have to manually delete the temporary directory ...")
        # sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
    return outfilename

def _write_script(scriptfilename,
                  wmin,wmax,dw,
                  costheta,
                  modelfilename,
                  marcsfile,
                  modelopacname,
                  metals,
                  alphafe,
                  indiv_abu, # dictionary with atomic number, abundance
                  vmicro,
                  spherical,
                  resultfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=False):
    """Write the script file for babsma and bsyn"""
    with open(scriptfilename,'w') as scriptfile:
        if bsyn:
            scriptfile.write("'NLTE : '  '.false.'\n")
            scriptfile.write("'NLTEINFOFILE : '  '{}/SPECIES_LTE_NLTE_00000000.dat'\n".format(TSDEPCOEFF_PATH))
        scriptfile.write("'LAMBDA_MIN:'  '%.3f'\n" % wmin)
        scriptfile.write("'LAMBDA_MAX:'  '%.3f'\n" % wmax)
        scriptfile.write("'LAMBDA_STEP:' '%.3f'\n" % dw)
        if bsyn:
            scriptfile.write("'INTENSITY/FLUX:' 'Flux'\n")
            scriptfile.write("'COS(THETA)    :' '%.3f'\n" % costheta)
            scriptfile.write("'ABFIND        :' '.false.'\n")
        if not bsyn:
            scriptfile.write("'MODELINPUT:' '%s'\n" % modelfilename)
        if marcsfile:
            scriptfile.write("'MARCS-FILE:' '.true.'\n")
        else:
            scriptfile.write("'MARCS-FILE:' '.false.'\n")
        scriptfile.write("'MODELOPAC:' '%s'\n" % modelopacname)
        if bsyn:
            scriptfile.write("'RESULTFILE :' '%s'\n"
                             % resultfilename)
        scriptfile.write("'METALLICITY:'    '%.3f'\n" % metals)
        scriptfile.write("'ALPHA/Fe   :'    '%.3f'\n" % alphafe)
        scriptfile.write("'HELIUM     :'    '0.00'\n")
        scriptfile.write("'R-PROCESS  :'    '0.00'\n")
        scriptfile.write("'S-PROCESS  :'    '0.00'\n")
        # Individual abundances: specify everything
        # It matters for continuum opacities too
        abundances = {}
        for Z, abund in solar_abundances_Z.items():
            if Z in [1,2]: continue # skip H and He
            ## TSFitPy does NOT do this
            #elif Z in [8, 10, 12, 14, 16, 18, 20, 22]:
            #    # According to Gustaffson+2008
            #    # The abundances of so-called α elements O, Ne, Mg, Si, S, Ar, Ca, and Ti, so I modify those
            #    abundances[Z] = abund + metals + alphafe
            else: abundances[Z] = abund + metals
        # Update using any specified individual abundances
        for Z, abund in indiv_abu.items():
            if Z in abundances: abundances[Z] = abund
            else: raise ValueError("indiv_abu has Z not in solar! {Z}")
        nabu= len(abundances)
        #if nabu > 0:
        scriptfile.write("'INDIVIDUAL ABUNDANCES:'   '%i'\n" % nabu)
        for abu in abundances:
            scriptfile.write("%i %.3f\n" % (abu,abundances[abu]))
        if bsyn:
            niso= len(isotopes)
            if niso > 0:
                scriptfile.write("'ISOTOPES : ' '%i'\n" % niso)
                for iso in isotopes:
                    scriptfile.write('%s %s\n' % (iso,isotopes[iso]))
            # Linelists
            nlines= len(linelistfilenames)
            scriptfile.write("'NFILES   :' '%i'\n" % nlines)
            for linelistfilename in linelistfilenames:
                scriptfile.write("%s\n" % linelistfilename)
            if spherical:
                scriptfile.write("'SPHERICAL:'  'T'\n")
            else:
                scriptfile.write("'SPHERICAL:'  'F'\n")
            scriptfile.write("30\n")
            scriptfile.write("300.00\n")
            scriptfile.write("15\n")
            scriptfile.write("1.30\n")
        # else:
        scriptfile.write("'XIFIX:' 'T'\n")
        scriptfile.write("%.3f\n" % vmicro)
    return None

def parse_model_atmosphere_file_params(model_atmosphere_file):
    with open(model_atmosphere_file,"r") as fp:
        pass
    header, _ = marcs.parse_marcs_model(model_atmosphere_file)
    Teff = header["Teff"]
    logg = header["logg"]
    vt = header["vturb"]
    MH = header["feh"]
    aFe = header["alphafe"]
    spherical = header["spherical"]
    return Teff, logg, vt, MH, aFe, spherical

def get_default_linelist_filenames(include_H=True, wmin=None, wmax=None):
    """
    Default linelists based on TSFitPy. This is GES for 4200-9200, VALD fills 3700-4200 and 9200-9800.
    I also added molecular lines for C to GES.
    https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/?p=%2Flinelist&mode=list
    
    TODO use wmin and wmax to cut lines appropriately
    TODO add a check for maximum number of lines allowable by Turbospectrum as it gets to be a lot
    """
    fnames = ["nlte_ges_linelist_jmg04sep2023_I_II",
        "vald-3700-3800-for-grid-nlte-04sep2023",
        #"vald-3800-4200-for-grid-nlte-04sep2023", # This one is broken! Not sure why.
        "vald-9200-9300-for-grid-nlte-04sep2023",
        "vald-9300-9800-for-grid-nlte-04sep2023",
        # "12C12C_GESv5.bsyn",
        # "12C13C_GESv5.bsyn",
        # "12C14N_GESv5.bsyn",
        # "12C15N_GESv5.bsyn",
        #"12CH_GESv5.bsyn",
        # "13C13C_GESv5.bsyn",
        # "13C14N_GESv5.bsyn",
        #"13CH_GESv5.bsyn",
        # "14NH_GESv5.bsyn",
    ]
    if include_H:
        fnames.append("Hlinedata")
    return [os.path.join(TSLINELIST_PATH, x) for x in fnames]

def get_vald_linelist_filenames(include_H=True):
    """
    Default VALD linelists from TSFitPy
    Searches in path TSLINELIST_PATH
    """
    fnames = ["vald-3700-3800-for-grid.list",
        "vald-3800-4300-for-grid.list",
        "vald-4300-4800-for-grid.list",
        "vald-4800-5300-for-grid.list",
        "vald-5300-5800-for-grid.list",
        "vald-5800-6300-for-grid.list",
        "vald-6300-6800-for-grid.list",
        "vald-6800-7300-for-grid.list",
        "vald-7300-7800-for-grid.list",
        "vald-7800-8300-for-grid.list",
        "vald-8300-8800-for-grid.list",
        "vald-8800-9300-for-grid.list",
        "vald-9300-9800-for-grid.list",
    ]
    if include_H:
        fnames.append("Hlinedata")
    return [os.path.join("/Users/alexji/lib/tssynth/data/linelists_vald/", x) for x in fnames]
