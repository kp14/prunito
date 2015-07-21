# -*- coding: utf-8 -*-
from setuptools import setup

def readme():
    with open('README') as f:
        return f.read()

setup(name='biocuration',
      version='0.6.2',
      description='A collection of tools for protein biocuration.',
      long_description=readme(),
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
      packages=['biocuration', 'biocuration/uniprotkb'],
      install_requires=['requests'],
      include_package_data=True)