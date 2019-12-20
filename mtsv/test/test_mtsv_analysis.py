import os
from mtsv.scripts.mtsv_analyze import *
import pytest
import datetime
import tables
import pandas as pd

@pytest.fixture(scope="function")
def existing_empty_datastore(tmpdir_factory):
    fn = tmpdir_factory.mktemp("datastores").join("empty_datastore.h5")
    store = pd.HDFStore(fn)
    store.close()
    return fn

@pytest.fixture(scope="function")
def existing_datastore(tmpdir_factory):
    fn = tmpdir_factory.mktemp("datastores").join("datastore.h5")
    ds = DataStore(fn)
    ds.add_dataset(
        "/datasets/t1", pd.DataFrame([1,2,3,4]),
        {'meta1': 1, 'meta2': 2})
    ds.add_dataset(
        "/datasets/t2", pd.DataFrame([5,6,7,8]),
        {'meta1': 3, 'meta2': 4})
    ds.close()
    return fn


@pytest.fixture(scope="function")
def empty_dir(tmpdir_factory):
    fn = tmpdir_factory.mktemp("empty_dir").join("datastore.h5")
    return fn


def add_df_assert_helper(name, meta, ds, df, existing=0):
    assert str(ds.store[name]) == str(df)
    assert ds.store.get_storer(name).attrs.metadata == meta
    assert len(ds.store.keys()) == (existing + 1)
    assert name in ds.store.keys()

def add_df_helper(fp, name='/datasets/t3'):
    ds = DataStore(fp)
    df = pd.DataFrame([9,10,11,12])
    meta = {'meta1': 5, 'meta2': 6}
    ds.add_dataset(name, df, meta)
    return name, meta, ds, df


def test_add_df_to_datastore_with_existing(existing_datastore):
    name, meta, ds, df = add_df_helper(existing_datastore)
    add_df_assert_helper(name, meta, ds, df, existing=2)

def test_add_df_to_empty_dir(empty_dir):
    name, meta, ds, df = add_df_helper(empty_dir)
    add_df_assert_helper(name, meta, ds, df)

def test_add_df_to_empty_datastore(existing_empty_datastore):
    name, meta, ds, df = add_df_helper(existing_empty_datastore)
    add_df_assert_helper(name, meta, ds, df)

def test_add_existing_df_to_existing_datastore_exception(existing_datastore):
    with pytest.raises(KeyError):
        add_df_helper(existing_datastore, name="/datasets/t1")
        
def test_get_dataset_from_existing_ds(existing_datastore):
    ds = DataStore(existing_datastore)
    result = ds.get_dataset("/datasets/t1")
    assert str(result) == str(pd.DataFrame([1, 2, 3, 4]))

def test_get_dataset_from_existing_ds_keyerror_exception(existing_datastore):
    ds = DataStore(existing_datastore)
    with pytest.raises(KeyError):
        ds.get_dataset("/datasets/t3")


def test_update_metadata_replace_keys(existing_datastore):
    ds = DataStore(existing_datastore)
    name = "/datasets/t1"
    meta = {'meta1': 1, 'meta2': 2}
    ds.update_metadata(name, meta)
    assert ds.store.get_storer(name).attrs.metadata == meta


def test_update_metadata_added_keys(existing_datastore):
    ds = DataStore(existing_datastore)
    name = "/datasets/t1"
    meta = {'meta3': 3}
    ds.update_metadata(name, meta)
    assert len(ds.store.get_storer(name).attrs.metadata) == 3

def test_get_metadata_from_data_store(existing_datastore):
    ds = DataStore(existing_datastore)
    name = "/datasets/t1"
    meta = ds.get_metadata(name)
    assert meta == {'meta1': 1, 'meta2': 2}

def test_key_error_raised_from_get_metadata(existing_empty_datastore):
    ds = DataStore(existing_empty_datastore)
    name = "/datasets/t1"
    with pytest.raises(KeyError):
        ds.get_metadata(name)

def test_key_error_raised_from_update_metadata(existing_empty_datastore):
    ds = DataStore(existing_empty_datastore)
    name = "/datasets/t1"
    with pytest.raises(KeyError):
        ds.update_metadata(name, {'meta1': 1})


def test_remove_dataset_from_datastore(existing_datastore):
    name = "/datasets/t1"
    ds = DataStore(existing_datastore)
    ds.remove(name)
    assert len(ds.store.keys()) == 1
    assert name not in ds.store.keys()

def test_remove_nonexistant_dataset_from_datastore_raises_keyerror(
    existing_empty_datastore):
    name = "/datasets/t1"
    ds = DataStore(existing_empty_datastore)
    with pytest.raises(KeyError):
        ds.remove(name)

def test_datastore_group_info(existing_datastore):
    ds = DataStore(existing_datastore)
    info = ds.group_info()
    expected = {
        "/datasets/t1": {
            "metadata": {'meta1': 1, 'meta2': 2},
            "size": ds.store['/datasets/t1'].memory_usage().sum()},
        "/datasets/t2": {
            "metadata": {'meta1': 3, 'meta2': 4},
            "size": ds.store['/datasets/t2'].memory_usage().sum()
        }}
    assert info == expected

def test_datastore_group_info_on_empty(existing_empty_datastore):
    ds = DataStore(existing_empty_datastore)
    info = ds.group_info()
    assert info == {}


def test_datastore_in_nonwritable_dir(existing_datastore):
    os.chmod(existing_datastore, 0o444)
    with pytest.raises(tables.FileModeError):
        add_df_helper(existing_datastore)


