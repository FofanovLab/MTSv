#!/usr/bin/env python
import os
from setuptools import setup, find_packages
from mtsv import VERSION

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


entry_points = '''
    [console_scripts]
    mtsv=mtsv.__main__:cli
    mtsv_setup=mtsv.mtsv_prep.main:main
'''

setup(name='mtsv',
      version=VERSION,
      description='Metagenomic analysis pipeline',
      author='Tara Furstenau',
      author_email='Tara.Furstenau@nau.edu',
      url='https://github.com/FofanovLab/MTSv',
      long_description=read('README.md'),
      install_requires=[
        'click',],
      packages=find_packages(),
      license='LICENSE.txt',
      entry_points=entry_points,
      include_package_data=True,
      package_data={
      'mtsv': [
          'commands/cmd_specs/*.yml',
          'scripts/*.py',
          'templates/*.html'
          'commands/snakefiles/*snek',
          'commands/snakefiles/Snakefile']
      },
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License']
      )



