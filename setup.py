#! python
"""Setup ms2rescore."""

from setuptools import setup


def get_version(path):
    """Get __version__ from Python file."""
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("__version__ = "):
                return line.strip().split(" = ")[1].strip("\"'")


setup(
    name="ms2rescore",
    version=get_version("ms2rescore/_version.py"),
    description="MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and retention times.",
    author="Ana Sílvia C. Silva, Ralf Gabriels, Tim Van Den Bossche",
    author_email="compomics.list@gmail.com",
    url="https://compomics.github.io/projects/ms2rescore/",
    packages=["ms2rescore"],
    include_package_data=True,
    entry_points={"console_scripts": ["ms2rescore=ms2rescore.__main__:main"],},
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    install_requires=[
        "importlib-resources;python_version<'3.7'",
        "numpy>=1.16.0,<2",
        "pandas>=0.24.0,<2",
        "scipy>=1.2.0,<2",
        "scikit-learn>=0.20.0,<1",
        "tqdm>=4.31.0,<5",
        "pyteomics>=4.1.0,<5",
        "xmltodict>=0.12.0,<1",
        "lxml>=4.5,<5",
        "ms2pip>=3.6,<4",
        "click>=7,<8",
        "deeplc>=0.1.17"
        "cascade-config",
    ],
    test_suite="ms2rescore.tests.test_ms2rescore.Tests",
    tests_require=["ms2rescore", "pytest>=4.3.0,<5",],
)
