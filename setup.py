import os
from setuptools import find_packages, setup


def extract_version():
    """
    Extracts version values from the main matplotlib __init__.py and
    returns them as a dictionary.
    """
    with open('nlabpy/__init__.py') as fd:
        for line in fd.readlines():
            if (line.startswith('__version__')):
                exec(line.strip())
    return locals()["__version__"]


def get_package_data():
    baseline_images = [
        'tests/baseline_images/%s/*' % x
        for x in os.listdir('tests/baseline_images')]

    return {
        'nlabpy':
        baseline_images +
        [
            "examples/*.html",
            "examples/*.txt"
        ]}

setup(
    name="nlabpy",
    version=extract_version(),
    author="Ilya Shamovsky",
    author_email="ilya.shamovsky@gmail.com",
    url="https://github.com/eco32i/nlabpy",
    license="MIT",
    packages=find_packages(),
    package_dir={"nlabpy": "nlabpy"},
    package_data=get_package_data(),
    description="scripts, functions and vizualization tools for NGS data",
    # run pandoc --from=markdown --to=rst --output=README.rst README.md
    long_description=open("README.rst").read(),
    # numpy is here to make installing easier... Needs to be at the
    # last position, as that's the first installed with
    # "python setup.py install"
    install_requires=["bokeh",
                      "pysam",
                      "pandas >= 0.16.0",
                      "numpy"],
    classifiers=['Intended Audience :: Science/Research',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Visualization',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5'],
    zip_safe=False)
