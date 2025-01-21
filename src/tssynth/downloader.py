import requests
from tqdm import tqdm
import os
import yaml
from importlib.resources import files as resource_files
import zipfile
import time

TSDEPCOEFF_PATH = os.environ.get('TSDEPCOEFF_PATH', None)
ALLMARCS_PATH = os.environ.get('ALLMARCS_PATH', None)

def download_file(url, local_path):
    """
    Downloads a file from the given URL and saves it to the specified local path.

    :param url: URL of the file to download
    :param local_path: Path where the downloaded file will be saved
    """
    response = requests.get(url, allow_redirects=True, stream=True)
    response.raise_for_status()
    
    total_size = int(response.headers.get('content-length', 0))
    block_size = 8192
    
    with open(local_path, 'wb') as file, tqdm(
        total=total_size, unit='iB', unit_scale=True
    ) as progress_bar:
        for chunk in response.iter_content(chunk_size=block_size):
            file.write(chunk)
            progress_bar.update(len(chunk))


def get_nlte_depgrid_info():
    ## TODO: this is not recommended, but let's keep it for now...
    nlte_info_path = resource_files('tssynth').joinpath('../../data/nlte_info.yml')
    with open(nlte_info_path, "r") as fp:
        nlte_info = yaml.load(fp, Loader=yaml.FullLoader)
    return nlte_info

def get_marcs_model_list():
    ## TODO: this is not recommended, but let's keep it for now...
    marcs_info_path = resource_files('tssynth').joinpath('../../data/model_list.txt')
    with open(marcs_info_path, "r") as fp:
        lines = [x.strip() for x in fp.readlines()]
    return lines

def download_nlte_depgrid(element):
    """
    Downloads the NLTE departure coefficients for the given element from the MPIA Bergemann Group
    From their Keeper.

    :param element: Element for which to download the NLTE departure coefficients
            See the list of available elements in tssynth/data/nlte_info.yml
    """
    if TSDEPCOEFF_PATH is None:
        raise ValueError("Environment variable TSDEPCOEFF_PATH is not set.")
    nlte_info = get_nlte_depgrid_info()
    if element not in nlte_info:
        raise ValueError(f"Element {element} not found in NLTE info file:\n{list(nlte_info.keys())}")
    element_dir = os.path.join(TSDEPCOEFF_PATH, element)
    os.makedirs(element_dir, exist_ok=True)
    all_files = nlte_info[element]
    for file in all_files:
        local_path = os.path.join(TSDEPCOEFF_PATH, os.path.join(element,file))
        if os.path.exists(local_path):
            print(f"File {local_path} already exists. Skipping.")
            continue
        url = f"https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/files/?p=%2Fdep-grids%2F{element}%2F{file}&dl=1"
        print(f"Downloading {file} from {url} to {local_path} (May be very large!)")
        download_file(url, local_path)
        print(f"Downloaded {file} to {local_path}")

def download_model_atmospheres():
    url = "https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/files/?p=%2Fatmospheres%2Fmarcs_standard_comp.zip&dl=1"
    model_list = get_marcs_model_list()
    missing_files = [model for model in model_list if not os.path.exists(os.path.join(ALLMARCS_PATH, model))]
    if not missing_files:
        print(f"All files from model_list.txt are present in {ALLMARCS_PATH}.")
        print("No need to download files.")
        return
    print(f"{ALLMARCS_PATH} is missing {len(missing_files)} model atmospheres out of {len(model_list)}")
    print(f"Downloading the entire model atmosphere archive from {url} to {ALLMARCS_PATH}")
    
    local_path = os.path.join(ALLMARCS_PATH, "marcs_standard_comp.zip")
    download_file(url, local_path)
    print(f"Downloaded model atmospheres to {local_path}, now unzipping...")
    start_time = time.time()
    with zipfile.ZipFile(local_path, 'r') as zip_ref:
        zip_ref.extractall(ALLMARCS_PATH)
    end_time = time.time()
    print(f"Unzipping completed in {end_time - start_time:.2f} seconds.")