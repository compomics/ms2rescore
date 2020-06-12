#! python

from setuptools import setup

setup(
    name='ms2rescore',
    version='0.2.1',
    description='MS²ReScore: PSM rescoring with predicted MS² peak intensities.',
    author='Ana Sílvia C. Silva',
    author_email='anascsilva@vib-ugent.be',
    url='https://www.github.com/compomics/ms2rescore',
    packages=['ms2rescore'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['ms2rescore=ms2rescore:main'],
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
        'pandas>=0.24.0,<1',
        'scipy>=1.2.0,<2',
        'scikit-learn>=0.20.0,<1',
        'tqdm>=4.31.0,<5',
        'pyteomics>=4.1.0,<5',
        'xmltodict>=0.12.0,<1',
    ],
    test_suite='ms2rescore.tests.test_ms2rescore.Tests',
    tests_require=[
        'numpy>=1.16.0,<2',
        'pandas>=0.24.0,<1',
        'scipy>=1.2.0,<2',
        'scikit-learn>=0.20.0,<1',
        'tqdm>=4.31.0,<5',
        'pyteomics>=4.1.0,<5',
        'xmltodict>=0.12.0,<1',
        'pytest>=4.3.0,<5',
    ],
)
