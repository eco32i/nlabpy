nlabpy
======

Genomic data parsing and visualization scripts.

Project Overview
----------------

**nlabpy** is a Python package designed for parsing, processing, and visualizing genomic data, specifically tailored for Next-Generation Sequencing (NGS) workflows. It provides a collection of scripts and functions to handle common tasks such as demultiplexing, filtering, and generating coverage plots or heatmaps.

**Key Technologies:**

*   **Language:** Python (Compatible with 2.7+ and 3.4+)
*   **Data Handling:** ``pandas``, ``numpy``, ``pysam``
*   **Visualization:** ``bokeh``, ``plotnine``
*   **Documentation:** Sphinx

Installation
------------

The project is set up to be installed as an editable package.

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/eco32i/nlabpy

    # Install dependencies and the package in editable mode
    pip install -e .

Usage
-----

Python API
^^^^^^^^^^

Import modules directly in your scripts or Jupyter notebooks:

.. code-block:: python

    from nlabpy.parse import seq
    from nlabpy.plotting import heatmap

CLI Scripts
^^^^^^^^^^^

Some functionalities are exposed as scripts in ``nlabpy/bin``. For example, to run the demultiplexer:

.. code-block:: bash

    python nlabpy/bin/demux.py --help

See example notebooks in the ``examples`` directory for more detailed usage.

Structure
---------

*   **nlabpy/**: The main source directory.

    *   **bin/**: Contains executable scripts, such as ``demux.py``.
    *   **parse/**: Modules for parsing different data formats (e.g., ``seq.py``, ``demux.py``, ``fqfilter.py``, ``anno.py``).
    *   **plotting/**: Visualization tools (e.g., ``heatmap.py``, ``genomecov.py``, ``deeptools.py``).
    *   **utils/**: Utility helpers including ``mprun.py`` for multiprocessing support.

*   **examples/**: Jupyter notebooks demonstrating usage.
*   **docs/**: Documentation source files using Sphinx.

Documentation
-------------

To build the documentation (assuming Sphinx is installed):

.. code-block:: bash

    cd docs
    make html