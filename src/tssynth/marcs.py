import numpy as np
import os, re, time, glob
from astropy.table import Table

def compress_marcs_standard_models(output_directory):
    """
    Compresses the standard MARCS models into a single file.
    """
    N_layer = 56

    ALLMARCS_PATH = os.environ.get("ALLMARCS_PATH")
    fnames = np.sort(glob.glob(f"{ALLMARCS_PATH}/*_st_*.mod"))
    assert os.path.exists(output_directory), output_directory
    print(f"Compressing {len(fnames)} Standard MARCS models in {ALLMARCS_PATH} to {os.path.abspath(output_directory)}.")
    all_headers_spherical = []
    all_stellarparams_spherical = []
    all_model_structures_spherical = []
    all_headers_planeparallel = []
    all_stellarparams_planeparallel = []
    all_model_structures_planeparallel = []
    start = time.time()
    for i, fname in enumerate(fnames):
        try:
            header, model_structure = parse_marcs_model(fname)
        except Exception as e:
            print(f"======={i}=======")
            print("FAILED ON", fname)
            print(e)
        assert len(model_structure) == N_layer, len(model_structure)
        if header["spherical"]:
            assert os.path.basename(fname).startswith("s"), fname
            all_headers_spherical.append(list(header.values()))
            all_stellarparams_spherical.append([header["Teff"], header["logg"], header["feh"]])
            all_model_structures_spherical.append(model_structure)
        else:
            assert os.path.basename(fname).startswith("p"), fname
            all_headers_planeparallel.append(list(header.values()))
            all_stellarparams_planeparallel.append([header["Teff"], header["logg"], header["feh"]])
            all_model_structures_planeparallel.append(model_structure)
        if i % 1000 == 0 and i > 0:
            elapsed = time.time() - start
            print(f"Processed {i} models in {elapsed:.2f} seconds.")
        # break

    colnames_headers = list(header.keys())

    headers_spherical = Table(rows=all_headers_spherical, names=colnames_headers)
    stellarparams_spherical = np.array(all_stellarparams_spherical, dtype=float)
    N_spherical = len(all_model_structures_spherical)
    model_structures_spherical = np.zeros((N_spherical, N_layer), dtype=all_model_structures_spherical[0].dtype)
    for i, model_structure in enumerate(all_model_structures_spherical):
        model_structures_spherical[i] = model_structure

    headers_planeparallel = Table(rows=all_headers_planeparallel, names=colnames_headers)
    stellarparams_planeparallel = np.array(all_stellarparams_planeparallel, dtype=float)
    N_planeparallel = len(all_model_structures_planeparallel)
    model_structures_planeparallel = np.zeros((N_planeparallel, N_layer), dtype=all_model_structures_planeparallel[0].dtype)
    for i, model_structure in enumerate(all_model_structures_planeparallel):
        model_structures_planeparallel[i] = model_structure

    np.savez_compressed(
        os.path.join(output_directory, "marcs_standard_models.npz"),
        headers_spherical=headers_spherical,
        model_structures_spherical=model_structures_spherical,
        stellarparams_spherical=stellarparams_spherical,
        headers_planeparallel=headers_planeparallel,
        model_structures_planeparallel=model_structures_planeparallel,
        stellarparams_planeparallel=stellarparams_planeparallel
    )

    
def parse_marcs_model(fname, get_all=False):
    with open(fname, "r") as fp:
        lines = [x[:-1] for x in fp.readlines()] # removes \n
    header = {
        "model_filename": lines[0],
        "Teff": float(lines[1].split()[0]),
        "Flux": float(lines[2].split()[0]),
        "logg": round(np.log10(float(lines[3].split()[0])),2),
        "vturb": float(lines[4].split()[0]),
        "mass": float(lines[5].split()[0]),
        "feh": float(lines[6].split()[0]),
        "alphafe": float(lines[6].split()[1]),
        "radius": float(lines[7].split()[0]),
        "luminosity": float(lines[8].split()[0]),
        "conv_alpha": float(lines[9].split()[0]),
        "conv_nu": float(lines[9].split()[1]),
        "conv_y": float(lines[9].split()[2]),
        "conv_beta": float(lines[9].split()[3]),
        "X": float(lines[10].split()[0]),
        "Y": float(lines[10].split()[1]),
        "Z": float(lines[10].split()[2]),
        "C13C12": float(re.search(r'12C/13C=(\d+)', lines[10]).group(1))
    }
    header["spherical"] = header["radius"] > 1
    assert lines[11].startswith("Logarithmic chemical number abundances, H always 12.00")
    abundances = {}
    for i, line in enumerate(lines[12:12 + 10]):
        values = list(map(float, line.split()))
        for j, value in enumerate(values):
            abundances[i * 10 + j + 1] = value
    assert len(abundances) == 92, len(abundances)

    depth_points_index = next(i for i, line in enumerate(lines) if re.match(r'^\d+ Number of depth points$', line.strip()))
    assert depth_points_index == 12 + 10, depth_points_index
    num_depth_points = int(lines[depth_points_index].split()[0])
    assert lines[depth_points_index + 1].strip() == "Model structure"
    cols_1_index = depth_points_index + 2
    cols_2_index = depth_points_index + 2 + num_depth_points + 1
    cols_1 = lines[cols_1_index].split()
    assert len(cols_1) == 9, cols_1
    cols_2 = lines[cols_2_index].split()
    assert cols_2[0] == "k"
    cols_2[0] = "k2"
    assert cols_2[1] == "lgTauR"
    cols_2[1] = "lgTauR2"
    assert len(cols_2) == 8, cols_2

    # Extracting the model structure
    def getfloat(x):
        try: return float(x)
        except Exception as e:
            if x == "******": return np.nan # this happens sometimes for mu
            else: raise e
    def formatted_read_1(line):
        try:
            widths = [3, 6, 8, 11, 8, 11, 11, 11, 11]
            return [getfloat(line[sum(widths[:i]):sum(widths[:i+1])].strip()) for i in range(len(widths))]
        except Exception as e: # There are a handful of MARCS models that don't follow the right widths
            return list(map(float, line.split()))
    def formatted_read_2(line):
        try:
            widths = [3, 6, 11, 11, 6, 11, 8, 14]
            return [getfloat(line[sum(widths[:i]):sum(widths[:i+1])].strip()) for i in range(len(widths))]
        except Exception as e: # There are a handful of MARCS models that don't follow the right widths
            return list(map(float, line.split()))

    model_structure_1 = np.array([formatted_read_1(line) for line in lines[cols_1_index + 1:cols_1_index + 1 + num_depth_points]])
    model_structure_2 = np.array([formatted_read_2(line) for line in lines[cols_2_index + 1:cols_2_index + 1 + num_depth_points]])
    
    # Concatenating model_structure_1 and model_structure_2 into a structured array
    dtype = [(name, float) for name in cols_1 + cols_2]
    model_structure = np.zeros(num_depth_points, dtype=dtype)
    for i, name in enumerate(cols_1):
        model_structure[name] = model_structure_1[:, i]
    for i, name in enumerate(cols_2):
        model_structure[name] = model_structure_2[:, i]
    
    assert lines[cols_2_index + num_depth_points + 1].startswith("Assorted logarithmic partial pressures")
    partial_pressures_lines = lines[cols_2_index + num_depth_points + 1:]
    
    if get_all:
        return header, model_structure, abundances, partial_pressures_lines
    else:
        return header, model_structure
    

def interpolate_marcs_model(Teff, logg, MH):
    """
    Interpolates from the MARCS model atmosphere grid to a new set of parameters.
    Only works for "st" (standard) MARCS models, and assumes vt = 2.0.
    
    This calls the Fortran program `interpol_modeles_nlte` which is in Turbospectrum_NLTE/interpolator.
    Requires the environment variable MARCS_GRID to be set to the directory containing the MARCS grid.
    """
    MARCS_GRID = os.environ.get("MARCS_GRID")

    raise NotImplementedError("This function is not implemented yet.")
