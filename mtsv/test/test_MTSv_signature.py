import pytest
from mtsv.scripts.MTSv_signature import *

def id_generator(start=1):
    i = start
    while True:
        yield i
        i += 1

ID_GEN = id_generator(1)

TEST_DATA_TUPLE = (
    ((294,300,303,28151,32002,36873,72557,76759,80842,85698,86265,158899,198620,200450,237609,262324,279113,587753,748247,930166,1259844,1274359,1283291,1458425,1458426,1534110,1662285,1785145,1842540
    ), "family", None),
    ((294,300,303,28151,32002,36873,72557,76759,80842,85698,86265,158899,198620,200450,237609,262324,279113,587753,748247,930166,1259844,1274359,1283291,1458425,1458426,1534110,1662285,1785145,1842540
    ), "genus", None),
    ((1392,1396,1428,658666), "genus", 1386),
    ((1392,1396,1428,658666), "family", 186817),
    ((1392,1396,1428,155322,580165,658666,1868655,1892404), "genus", 1386),
    ((1392,1396,1428,155322,580165,658666,1868655,1892404), "family", 186817),
    ((1392,1396,1428,658666,240303), "genus", None),
    ((1392,1396,1428,658666,240303), "family", 186817))

TEST_DATA_LIST = [
    [294,300,303,28151,32002,36873,72557,76759,80842,85698,86265,158899,198620,200450,237609,262324,279113,587753,748247,930166,1259844,1274359,1283291,1458425,1458426,1534110,1662285,1785145,1842540],
    [248026],
    [1392,1396,1428,658666],
    [1392,1396,1428,155322,580165,658666,1868655,1892404],
    [34018],
    [1392,1396,1428,658666,240303],
    [1392,1396,1428,658666]
]



def format_sig_line(taxids, counts="1_1_1"):
    return "R{0}_{1}:{2}\n".format(
        next(ID_GEN), counts, ",".join([str(tax) for tax in taxids]))

@pytest.fixture(scope="module")
def merge_file(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("merged.clp")
    out_string = ""
    for taxids in TEST_DATA_LIST:
        out_string += format_sig_line(taxids)
    fn.write(out_string)
    return fn
        
@pytest.mark.parametrize("taxids, rank, expected", TEST_DATA_TUPLE)
def test_get_common_ancestor(taxids, rank, expected):
    assert get_common_ancestor(taxids, rank) == expected

def test_separate_multiples_from_singletons(merge_file):
    result = separate_multiples_from_singletons(merge_file)
    for singleton, expected in zip(result[0], [b"R2", b"R5"]):
        assert singleton.read_id == expected
    for e, (k, v) in zip(
        [[b"R1"], [b"R3", b"R7"], [b"R4"], [b"R6"]], result[1].items()):
        for record, exp in zip(v, e):
            assert record.read_id == exp
        
    