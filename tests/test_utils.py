from tssynth import utils
import os, shutil

def test_element_to_atomic_number():
    assert utils.element_to_atomic_number("H") == 1, utils.element_to_atomic_number("H")
    assert utils.element_to_atomic_number("Fe") == 26
    assert utils.element_to_atomic_number("U") == 92
    assert utils.element_to_atomic_number("X") == None

def test_parse_XFe_dict():
    assert utils.parse_XFe_dict({"H": 0.0, "He": 0.0, "C": 0.0}) == {1: 0.0, 2: 0.0, 6: 0.0}
    assert utils.parse_XFe_dict({1: 0.0, 2: 0.0, 6: 0.0}) == {1: 0.0, 2: 0.0, 6: 0.0}

def test_mkdtemp():
    twd = utils.mkdtemp()
    assert os.path.exists(twd)
    shutil.rmtree(twd)
