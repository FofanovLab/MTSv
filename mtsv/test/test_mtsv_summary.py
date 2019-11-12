from ete3 import NCBITaxa
import pytest
import os
import pandas as pd
from mtsv.scripts.mtsv_summary import *

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

MERGED_DATA = os.path.join(TEST_DIR, 'merged.clp')

TEST_DATA_LIST = [
    [294, 300, 303, 28151, 32002, 36873, 72557, 76759, 80842, 85698, 86265, 158899, 198620, 200450, 237609, 262324,
        279113, 587753, 748247, 930166, 1259844, 1274359, 1283291, 1458425, 1458426, 1534110, 1662285, 1785145, 1842540],
    [248026],
    [1392, 1396, 1428, 658666],
    [1392, 1396, 1428, 155322, 580165, 658666, 1868655, 1892404],
    [34018],
    [1392, 1396, 1428, 658666, 240303],
    [1392, 1396, 1428, 658666]
]


lineages = [
    [1392, 1386, 186817, 1385, 91061, 1239, 2],
    [222, 222, 506, 80840, 28216, 1224, 2],
    [32630, 32630, 32630, 32630, 32630, 32630, 32630]]

taxa = [l[0] for l in lineages]


lineage_table_lineage_level_test_data = list(zip(taxa, lineages))

lineage_level_data = (
    (taxa, level,
    list(set([lineage[i] for lineage in lineages])))
    for i,level in enumerate(LEVELS)) 

lineage_expanded_data = [[1392], [1392,222], [32630], [1392, 222, 32630]]



def id_generator(start=1):
    i = start
    while True:
        yield i
        i += 1


ID_GEN = id_generator(1)


def format_sig_line(taxids, counts="10"):
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


@pytest.fixture(scope="module")
def data_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("data")


@pytest.fixture(scope="module")
def ncbi():
    ncbi = NCBITaxa()
    yield ncbi


@pytest.fixture(scope="module")
def lineage_table():
    df = pd.DataFrame(lineages,
        columns = LEVELS, index=taxa)
    yield df


@pytest.fixture(scope="module")
def lineage_instance(ncbi):
    lineage = LineageTable(ncbi)
    df = pd.DataFrame(lineages,
                      columns=LEVELS, index=taxa)
    lineage.lineage_table = df
    yield lineage

@pytest.mark.parametrize("taxon, expected", lineage_table_lineage_level_test_data)
def test_LineageTable_get_lineage_levels(taxon, expected, ncbi):
    lineage = LineageTable(ncbi)
    assert lineage.get_lineage_levels(taxon) == expected


def test_LineageTable_update_lineage_table_does_not_add_existing_taxa(
        lineage_instance):
    lineage_instance.lineage_table = lineage_table    
    shape1 = lineage_instance.lineage_table.shape
    lineage_instance.update_lineage_table(pd.Series([[1392]]))
    assert shape1 == lineage_instance.lineage_table.shape
    
def test_LineageTable_update_lineage_table_does_add_taxa(ncbi):
    lineage = LineageTable(ncbi)
    lineage.update_lineage_table(TEST_DATA_LIST)
    unique = set()
    for l in TEST_DATA_LIST:
        unique.update(l)
    assert list(unique) == list(lineage.lineage_table.index)


@pytest.mark.parametrize("taxa, level, expected", lineage_level_data)
def test_LineageTable_get_level_from_taxa(
        taxa, level, expected, lineage_instance):
    lineage_instance.lineage_table=lineage_table
    result = lineage_instance.get_level_from_taxa(taxa, level)
    assert list(set(result)) == expected


def test_LineageTable__get_new_unique_taxa(lineage_instance):
    lineage_table_shape1 = lineage_instance.lineage_table.shape
    lineage_table._get_new_unique_taxa(
        [[l] for l in lineage_instance.lineage_table.index])
    assert lineage_instance.lineage_table.shape == lineage_table_shape1


@pytest.mark.parametrize("taxa, level, expected", lineage_expanded_data)
def test_expand_taxa_array(ncbi):
    lineage = LineageTable(ncbi)
    lineage.update_lineage_table(TEST_DATA_LIST)
    result = lineage.expand_taxa_array(taxa, level)
    assert result == expected
