import configparser
from mtsv.commands import (
    init, readprep, binning, summary, analyze, extract, pipeline)
from mtsv.parsing import specfile_read
import argutils
 #from mtsv import specfile


_commands = [
    'readprep',
    'binning',
    'analyze',
    'summary',
    'extract'
]


def create_config_file():
    cfg_file_str = ""
    for cmd in _commands:
        spec = specfile_read(cmd)
        opts = argutils.read.from_yaml(spec)
        if opts:
            cfg_file_str += argutils.export.to_config(
                cmd.upper(), opts) + "\n"
    return cfg_file_str



Init = init.Init
Readprep = readprep.Readprep
Binning = binning.Binning
Summary = summary.Summary
Analyze = analyze.Analyze
Extract = extract.Extract
Pipeline = pipeline.Pipeline

