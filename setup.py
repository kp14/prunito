# -*- coding: utf-8 -*-
import codecs
import os
import re

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    # Taken from pip setup.py
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    # Taken from pip setup.py
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

long_description = read('README')

setup(name='biocuration',
      version=find_version('biocuration', '__init__.py'),
      description='A collection of tools for protein biocuration.',
      long_description=long_description,
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: Public Domain',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization',
      ],
      url='https://bitbucket.org/kp14/biocuration-project',
      author='kp14',
      author_email='',
      license='Public Domain',
      packages=['biocuration', 'biocuration/uniprot', 'biocuration/distil'],
      install_requires=['requests',
                        'lxml',
                        'beautifulsoup4',
                        'venndy',
                        ],
      include_package_data=True)