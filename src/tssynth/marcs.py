import numpy as np
import os, re, time, glob
import subprocess
from astropy.table import Table
from . import utils

def compress_marcs_standard_models(output_directory):
    """
    Compresses the standard MARCS models into a single file.
    This is not used right now, but someday I would like to replace
    the fortran interpolators with python interpolators, and this will
    be helpful for that.
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

def write_marcs_model(header, model_structure):
    raise NotImplementedError("This function is not implemented yet.")

def interpolate_marcs_model(Teff, logg, MH, outpath, spherical=True):
    """
    Runs the fortran interpolator (in TSINTERP_PATH) to interpolate the MARCS models.
    Looks in the ALLMARCS_PATH for the MARCS models.

    Parameters:
    -----------
    Teff : float
        Effective temperature of the star in Kelvin.
    logg : float
        Surface gravity of the star
    MH : float
        Metallicity [M/H] of the star.
    outpath : str
        Path where the interpolated model will be saved.
    spherical : bool, optional
        If True, use spherical MARCS models. If False, use plane-parallel MARCS models. Default is True.
        logg > 3.5 is not available for spherical models and logg < 3.0 is not available for plane-parallel models.
        We do not automatically specify which one, the user must choose correctly.
    Raises:
    -------
    ValueError
        If the input parameters are out of the valid range for MARCS models.
    Notes:
    ------
    This function calls the Fortran program `interpol_modeles` which is in Turbospectrum_NLTE/interpolator.
    Example:
    --------
    >>> interpolate_marcs_model(5777, 4.44, 0.0, "/path/to/output.interpol")
    """

    ALLMARCS_PATH = os.environ.get("ALLMARCS_PATH")

    ## Validate input parameters
    if spherical and logg > 3.5: raise ValueError("No spherical MARCS models for logg > 3.5.")
    if not spherical and logg < 3.0: raise ValueError("No plane-parallel MARCS models for logg < 3.0.")
    if logg > 5.5: raise ValueError("No MARCS models for logg > 5.5.")
    if logg < -0.5: raise ValueError("No MARCS models for logg < -0.5.")
    if MH > 0.5: raise ValueError("No MARCS models for MH > 0.5.")
    if MH < -5.0: raise ValueError("No MARCS models for MH < -5.0.")
    if Teff > 8000: raise ValueError("No MARCS models for Teff > 8000.")
    if Teff < 2500: raise ValueError("No MARCS models for Teff < 2500.")

    sstr = "s" if spherical else "p"
    massstr = "1.0" if spherical else "0.0"
    
    ## Search through input model atmospheres for models to interpolate
    fnames = np.sort(glob.glob(f"{ALLMARCS_PATH}/{sstr}*_t02_st_*.mod"))
    def parse_marcs_filenames(fname):
        """
        Parses the MARCS filenames to extract the stellar parameters.
        """
        parts = os.path.basename(fname).split("_")
        Teff = int(parts[0][1:])
        logg = float(parts[1][1:])
        MH = float(parts[5][1:])
        return Teff, logg, MH
    marcspoints = np.array([parse_marcs_filenames(fname) for fname in fnames])
    
    points = _find_surrounding_points(marcspoints, Teff, logg, MH)
    selected_files = []
    for point in points:
        Teff_point, logg_point, MH_point = point
        fname_start = f"{sstr}{Teff_point:4.0f}_g{logg_point:+4.1f}_m{massstr}_t02_st_z{MH_point:+5.2f}"
        for fname in fnames:
            if os.path.basename(fname).startswith(fname_start):
                selected_files.append(fname)
                break
        else:
            raise ValueError(f"Could not find {fname_start} in {ALLMARCS_PATH}.\n{points}")

    ## Run the fortran interpolator
    _run_interpolator_lte(Teff, logg, MH, selected_files,
                          outpath, verbose=False)
    
    assert os.path.exists(outpath), f"Interpolator did not create the output file {outpath}."
    return outpath

def _find_surrounding_points(marcspoints, Teff, logg, MH, max_expansions=5):
    """
    Find the 8 surrounding points for interpolation in a 3D grid.
    Needs to be a cuboid. Does a recursive search through the grid to find the smallest number of expansions that contains the points.
    This is slow if it needs to expand many times due to the recursion.

    Parameters:
    -----------
    marcspoints : ndarray
        Array of available grid points with shape (N, 3), where each row is [Teff, logg, MH].
    Teff : float
        Target effective temperature.
    logg : float
        Target surface gravity.
    MH : float
        Target metallicity.
    max_expansions : int, optional
        Maximum number of expansions to search for surrounding points. Default is 99.

    Returns:
    --------
    ndarray
        Array of 8 surrounding points, each row is [Teff, logg, MH].
    
    Raises:
    -------
    RuntimeError
        If 8 surrounding points are not found after the maximum number of expansions (default 3).
    """
    unique_Teffs = np.unique(marcspoints[:, 0])
    unique_loggs = np.unique(marcspoints[:, 1])
    unique_MHs = np.unique(marcspoints[:, 2])

    ## Check if the target point is within the grid boundaries
    if Teff < unique_Teffs[0] or Teff > unique_Teffs[-1]:
        raise ValueError(f"Teff={Teff} is outside the range of available MARCS models: {unique_Teffs[0]} to {unique_Teffs[-1]}.")
    if logg < unique_loggs[0] or logg > unique_loggs[-1]:
        raise ValueError(f"logg={logg} is outside the range of available MARCS models: {unique_loggs[0]} to {unique_loggs[-1]}.")
    if MH < unique_MHs[0] or MH > unique_MHs[-1]:
        raise ValueError(f"MH={MH} is outside the range of available MARCS models: {unique_MHs[0]} to {unique_MHs[-1]}.")

    def find_bounds(values, target):
        """Find lower and upper bounds for a target value in sorted values."""
        lower_idx = np.searchsorted(values, target, side="right") - 1
        upper_idx = lower_idx + 1

        lower_idx = max(lower_idx, 0)
        upper_idx = min(upper_idx, len(values) - 1)
        return lower_idx, upper_idx

    # Initial faces for Teff, logg, and MH
    Teff_lower_idx, Teff_upper_idx = find_bounds(unique_Teffs, Teff)
    logg_lower_idx, logg_upper_idx = find_bounds(unique_loggs, logg)
    MH_lower_idx, MH_upper_idx = find_bounds(unique_MHs, MH)
    indices = [Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx]

    ## Check that all the points are valid
    def check_indices_in_grid(indices):
        Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx = indices
        Teff_lower, Teff_upper = unique_Teffs[Teff_lower_idx], unique_Teffs[Teff_upper_idx]
        logg_lower, logg_upper = unique_loggs[logg_lower_idx], unique_loggs[logg_upper_idx]
        MH_lower, MH_upper = unique_MHs[MH_lower_idx], unique_MHs[MH_upper_idx]
        points = []
        for T in [Teff_lower, Teff_upper]:
            for g in [logg_lower, logg_upper]:
                for m in [MH_lower, MH_upper]:
                    iiT = marcspoints[:,0] == T
                    iig = marcspoints[:,1] == g
                    iiM = marcspoints[:,2] == m
                    ii = iiT & iig & iiM
                    if ii.sum() == 1: points.append((T, g, m))
        return len(points) == 8, np.array(points)
    valid, points = check_indices_in_grid(indices)
    if valid: return points

    ## Expand the cuboid until all points are found
    ## Do it recursively
    def find_best_cuboid(indices, N_left):
        # Finish cases: all points in grid or face out of bounds or max exceeded
        if N_left == 0: return 9999999, indices
        Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx = indices
        if Teff_lower_idx < 0 or Teff_upper_idx >= len(unique_Teffs) or \
           logg_lower_idx < 0 or logg_upper_idx >= len(unique_loggs) or \
           MH_lower_idx < 0 or MH_upper_idx >= len(unique_MHs):
            return 9999999, indices
        if check_indices_in_grid(indices)[0] == True: return 0, indices
        
        # Find the best expansion in all cases from this point
        N_Teff_lower, indices_Teff_lower = find_best_cuboid([Teff_lower_idx - 1, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx], N_left-1)
        N_Teff_upper, indices_Teff_upper = find_best_cuboid([Teff_lower_idx, Teff_upper_idx + 1, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx], N_left-1)
        N_logg_lower, indices_logg_lower = find_best_cuboid([Teff_lower_idx, Teff_upper_idx, logg_lower_idx - 1, logg_upper_idx, MH_lower_idx, MH_upper_idx], N_left-1)
        N_logg_upper, indices_logg_upper = find_best_cuboid([Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx + 1, MH_lower_idx, MH_upper_idx], N_left-1)
        N_MH_lower, indices_MH_lower = find_best_cuboid([Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx - 1, MH_upper_idx], N_left-1)
        N_MH_upper, indices_MH_upper = find_best_cuboid([Teff_lower_idx, Teff_upper_idx, logg_lower_idx, logg_upper_idx, MH_lower_idx, MH_upper_idx + 1], N_left-1)

        # These are sorted in order of preference, as we get the smallest index satisfying Ns==min(Ns)
        Ns = np.array([N_MH_upper, N_MH_lower, N_Teff_upper, N_Teff_lower, N_logg_upper, N_logg_lower])
        best_ix = np.where(np.min(Ns)==Ns)[0][0]
        if best_ix == 0: return N_MH_upper + 1, indices_MH_upper
        elif best_ix == 1: return N_MH_lower + 1, indices_MH_lower
        elif best_ix == 2: return N_Teff_upper + 1, indices_Teff_upper
        elif best_ix == 3: return N_Teff_lower + 1, indices_Teff_lower
        elif best_ix == 4: return N_logg_upper + 1, indices_logg_upper
        elif best_ix == 5: return N_logg_lower + 1, indices_logg_lower
        else: raise RuntimeError("This should not happen.")
    
    N_best, indices_best = find_best_cuboid(indices, max_expansions)
    print(f"Found best cuboid after {N_best} expansions.")
    if N_best >= 9999999: raise RuntimeError("Failed to find 8 surrounding points for {Teff}, {logg}, {MH}.")
    
    ## Triple check
    valid, points = check_indices_in_grid(indices_best)
    assert valid, points
    return points

def _run_interpolator_lte(Teff, logg, MH, marcs_model_list,
                        outpath, verbose=False):
    """
    Runs the Fortran interpolator to interpolate the MARCS models.
    Based on https://github.com/TSFitPy-developers/TSFitPy/blob/main/scripts/turbospectrum_class_nlte.py
            _interpolate_one_atmosphere
    """
    interp_exec_path = os.environ.get("TSINTERP_PATH")
    
    if verbose:
        stdout = None
        stderr = subprocess.STDOUT
    else:
        stdout = open('/dev/null', 'w')
        stderr = subprocess.STDOUT
    
    # Write configuration input for interpolator
    interpol_config = ""
    for marcs_model in marcs_model_list:
        assert os.path.exists(marcs_model), marcs_model
        interpol_config += "'{}'\n".format(marcs_model)
    interpol_config += "'{}'\n".format(outpath) # .interpol output file
    interpol_config += "'/dev/null'\n" # .alt output file, not needed
    interpol_config += "{}\n".format(Teff)
    interpol_config += "{}\n".format(logg)
    interpol_config += "{}\n".format(MH)
    interpol_config += ".false.\n"  
    interpol_config += ".false.\n"  
    interpol_config += "'/dev/null'\n" # .test output file, not needed

    # Now we run the FORTRAN model interpolator
    try:
        p = subprocess.Popen([os.path.join(interp_exec_path, 'interpol_modeles')],
                            stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
        p.stdin.write(bytes(interpol_config, 'utf-8'))
        stdout, stderr = p.communicate()
    except subprocess.CalledProcessError as e:
        print(e)
        raise RuntimeError("MARCS model atmosphere interpolation failed. Config:\n"+interpol_config)

    return outpath

def _run_interpolator_nlte(Teff, logg, MH, marcs_model_list):
    """
    Runs the Fortran interpolator for both the MARCS model atmospheres
    and the NLTE departure coefficient grids.
    Based on https://github.com/TSFitPy-developers/TSFitPy/blob/main/scripts/turbospectrum_class_nlte.py
            _interpolate_one_atmosphere
    """
    # for element in self.model_atom_file:
    #     element_abundance = self._get_element_abundance(element)
    #     # Write configuration input for interpolator
    #     interpol_config = ""
    #     for line in marcs_model_list:
    #         interpol_config += "'{}{}'\n".format(self.marcs_grid_path, line)
    #     interpol_config += "'{}.interpol'\n".format(output)
    #     interpol_config += "'{}.alt'\n".format(output)
    #     interpol_config += "'{}_{}_coef.dat'\n".format(output, element)  # needed for nlte interpolator
    #     interpol_config += "'{}'\n".format(os_path.join(self.departure_file_path, self.depart_bin_file[
    #         element]))  # needed for nlte interpolator
    #     interpol_config += "'{}'\n".format(os_path.join(self.departure_file_path, self.depart_aux_file[
    #         element]))  # needed for nlte interpolator
    #     interpol_config += "{}\n".format(self.aux_file_length_dict[element])
    #     interpol_config += "{}\n".format(self.t_eff)
    #     interpol_config += "{}\n".format(self.log_g)
    #     interpol_config += "{:.6f}\n".format(round(float(self.metallicity), 6))
    #     interpol_config += "{:.6f}\n".format(round(float(element_abundance), 6))
    #     interpol_config += ".false.\n"  # test option - set to .true. if you want to plot comparison model (model_test)
    #     interpol_config += ".false.\n"  # MARCS binary format (.true.) or MARCS ASCII web format (.false.)?
    #     interpol_config += "'{}.test'\n".format(output)

    #     # Now we run the FORTRAN model interpolator
    #     try:
    #         if self.atmosphere_dimension == "1D":
    #             p = subprocess.Popen([os_path.join(self.interpol_path, 'interpol_modeles_nlte')],
    #                                     stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
    #             p.stdin.write(bytes(interpol_config, 'utf-8'))
    #             stdout, stderr = p.communicate()
    #         elif self.atmosphere_dimension == "3D":
    #             p = subprocess.Popen([os_path.join(self.interpol_path, 'interpol_multi_nlte')],
    #                                     stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
    #             p.stdin.write(bytes(interpol_config, 'utf-8'))
    #             stdout, stderr = p.communicate()
    #     except subprocess.CalledProcessError:
    #         return {
    #             "interpol_config": interpol_config,
    #             "errors": "MARCS model atmosphere interpolation failed."
    #         }
    raise NotImplementedError("NLTE interpolation not implemented yet.")
    