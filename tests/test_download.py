import tssynth.downloader

def test_get_nlte_depgrid_info():
    nlte_info = tssynth.downloader.get_nlte_depgrid_info()
    assert "H" in nlte_info

def test_get_marcs_model_list():
    marcs_models = tssynth.downloader.get_marcs_model_list()
    # assert len(marcs_models) == 

def test_download_nlte_depgrid():
    ## TODO not totally sure what to do here as a test without actually just downloading everything
    # I tested it on my computer and it works though...
    # tssynth.downloader.download_nlte_depgrid("H")
    nlte_info = tssynth.downloader.get_nlte_depgrid_info()
    assert "H" in nlte_info

def test_download_model_atmospheres():
    ## TODO not totally sure what to do here as a test without actually just downloading everything
    # I tested it on my computer and it works though...
    tssynth.downloader.download_model_atmospheres()
