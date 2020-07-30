#! python
"""Setup ms2rescore."""

from setuptools import setup
from ms2rescore._version import __version__

setup(
    name='ms2rescore',
    version=__version__,
    description='MS²ReScore: PSM rescoring with predicted MS² peak intensities.',
    author='Ana Sílvia C. Silva, Ralf Gabriels, Tim Van Den Bossche',
    author_email='compomics.list@gmail.com',
    url='https://compomics.github.io/projects/ms2rescore/',
    packages=['ms2rescore'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['ms2rescore=ms2rescore.__main__:main'],
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta"
    ],
    install_requires=[
        'numpy>=1.16.0,<2',
        'pandas>=0.24.0,<2',
        'scipy>=1.2.0,<2',
        'scikit-learn>=0.20.0,<1',
        'tqdm>=4.31.0,<5',
        'pyteomics>=4.1.0,<5',
        'xmltodict>=0.12.0,<1',
        'lxml>=4.5,<5',
        'ms2pip>=3.6,<4',
        'click>=7,<8',
    ],
    test_suite='ms2rescore.tests.test_ms2rescore.Tests',
    tests_require=[
        'numpy>=1.16.0,<2',
        'pandas>=0.24.0,<2',
        'scipy>=1.2.0,<2',
        'scikit-learn>=0.20.0,<1',
        'tqdm>=4.31.0,<5',
        'pyteomics>=4.1.0,<5',
        'xmltodict>=0.12.0,<1',
        'lxml>=4.5,<5'
        'pytest>=4.3.0,<5',
    ],
)
