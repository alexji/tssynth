import os
import numpy as np
import re

def parse_marcs_model(fname, get_all=False):
    with open(fname, "r") as fp:
        lines = [x.strip() for x in fp.readlines()]
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
    header["spherical"] = header["mass"] > 0
    assert lines[11].startswith("Logarithmic chemical number abundances, H always 12.00")
    abundances = {}
    for i, line in enumerate(lines[12:12 + 10]):
        values = list(map(float, line.split()))
        for j, value in enumerate(values):
            abundances[i * 10 + j + 1] = value
    assert len(abundances) == 92, len(abundances)

    depth_points_index = next(i for i, line in enumerate(lines) if re.match(r'^\d+ Number of depth points$', line))
    assert depth_points_index == 12 + 10
    num_depth_points = int(lines[depth_points_index].split()[0])
    assert lines[depth_points_index + 1] == "Model structure"
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
    model_structure_1 = np.array([list(map(float, line.split())) for line in lines[cols_1_index + 1:cols_1_index + 1 + num_depth_points]])
    model_structure_2 = np.array([list(map(float, line.split())) for line in lines[cols_2_index + 1:cols_2_index + 1 + num_depth_points]])
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
