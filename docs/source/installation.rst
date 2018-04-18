.. _installation:

Installing prunito
==================

Prunito is packaged as a wheel and can be installed using `pip <https://pip.pypa.io/en/stable/>`_.
This will also install all the dependencies.
Python >=3.6 is required.
As usual, it probably best to install prunito and its dependencies into a dedicated virtual environment.

#. Create a virtual environment

    Using ``conda`` from the `Anaconda Python distribution <https://www.continuum.io/downloads>`_ :

        * Create a new environment (env) called *divvy* which runs Python 3.6 and has ``pip`` installed::

            conda create -n prunito python=3.6 pip

        * Activate the env::

            conda activate prunito

    Using *venv* in a regular Python installation. Python is usually available out of the box on Linux:

        * Ensure that the version of Python used is 3.6 or higher::

            python --version

        * Create a new environment (env) called *prunito*. ``ensurepip`` will bootstrap pip into the env::

            pyvenv /path/to/new/virtual/environment/prunito

#. Install prunito with its dependencies::

        pip install prunito-<version-py3-none-any>.whl

Dependencies
============

All of the following packages have to be installed some of which come with their own dependencies but those should
be taken care of by running the usual pip command:

* `requests <http://docs.python-requests.org/en/master/>`_
* `lxml <http://lxml.de/>`_
* `venndy <https://github.com/kp14/venndy>`_ (optional)
