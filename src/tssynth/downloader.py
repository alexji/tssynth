import requests
from tqdm import tqdm
def download_file(url, local_path):
    """
    Downloads a file from the given URL and saves it to the specified local path.

    :param url: URL of the file to download
    :param local_path: Path where the downloaded file will be saved
    """
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    total_size = int(response.headers.get('content-length', 0))
    block_size = 8192
    
    with open(local_path, 'wb') as file, tqdm(
        total=total_size, unit='iB', unit_scale=True
    ) as progress_bar:
        for chunk in response.iter_content(chunk_size=block_size):
            file.write(chunk)
            progress_bar.update(len(chunk))
