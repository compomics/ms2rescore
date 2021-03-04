#! python
"""Setup ms2rescore."""

from setuptools import setup


def get_version(path):
    """Get __version__ from Python file."""
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("__version__ = "):
                return line.strip().split(" = ")[1].strip("\"'")


def get_readme():
    with open("README.md", "r") as fh:
        readme = fh.read()
    return readme


setup(
    name="ms2rescore",
    version=get_version("ms2rescore/_version.py"),
    license="apache-2.0",
    description="MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and retention times.",
    long_description=get_readme(),
    long_description_content_type="text/markdown",
    author="Ana Sílvia C. Silva, Ralf Gabriels, Tim Van Den Bossche",
    author_email="compomics.list@gmail.com",
    url="https://compomics.github.io/projects/ms2rescore/",
    project_urls={
        "Documentation": "http://compomics.github.io/projects/ms2rescore",
        "Source": "https://github.com/compomics/ms2rescore",
        "Tracker": "https://github.com/compomics/ms2rescore/issues",
        "Publication": "https://doi.org/10.1093/bioinformatics/btz383",
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    keywords=[
        "MS2ReScore",
        "MS2PIP",
        "DeepLC",
        "Percolator",
        "Proteomics",
        "peptide",
        "peak intensity prediction",
        "spectrum",
        "machine learning",
    ],
    packages=["ms2rescore"],
    include_package_data=True,
    entry_points={"console_scripts": ["ms2rescore=ms2rescore.__main__:main"]},
    python_requires=">=3.7",
    install_requires=[
        "importlib-resources;python_version<'3.7'",
        "numpy>=1.16.0,<2",
        "pandas>=0.24.0,<2",
        "scikit-learn>=0.20.0,<1",
        "scipy>=1.2.0,<2",
        "tqdm>=4.31.0,<5",
        "pyteomics>=4.1.0,<5",
        "lxml>=4.5,<5",
        "ms2pip>=3.6.3,<4",
        "click>=7,<8",
        "cascade-config",
    ],
    extras_require={
        'deeplc': ["deeplc>=0.1.17"]
    },
    test_suite="tests",
    tests_require=["ms2rescore", "pytest>=4.3.0,<5"],
)
