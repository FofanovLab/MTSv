#!/usr/bin/env python
import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='mtsv',
      version='1.0',
      description='Metagenomic analysis pipeline',
      author='Tara Furstenau',
      author_email='Tara.Furstenau@nau.edu',
      url='https://github.com/FofanovLab/MTSv',
      long_description=read('README.md'),
      packages=find_packages(),
      license='LICENSE.txt',
      entry_points={'console_scripts': ['mtsv=mtsv.main:main',
      'mtsv_setup=mtsv.mtsv_prep.main:main',
      'mtsv_plugins=mtsv.mtsv_plugin.main:main']},
      include_package_data=True,
      package_data={
      'mtsv': [
          'commands/cmd_specs/*.yml',
          'scripts/*.py',
          'commands/snakefiles/*snek',
          'data/expected.data']
      },
      data_files=[
          ('mtsv/ext', ['mtsv/ext/target/release/mtsv-binner',
                   'mtsv/ext/target/release/mtsv-build',
                   'mtsv/ext/target/release/mtsv-chunk',
                   'mtsv/ext/target/release/mtsv-collapse',
                   'mtsv/ext/target/release/mtsv-readprep',
                   'mtsv/ext/target/release/mtsv-signature',
                   'mtsv/ext/target/release/mtsv-tree-build',
                   'mtsv/ext/mtsv-db-build'])
      ],
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License']
      )



