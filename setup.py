#! python

from setuptools import setup

setup(
    name='ms2rescore',
    version='0.2.1',
    description='MS2ReScore: PSM rescoring with predicted MS² peak intensities.',
    author='Ana Sílvia C. Silva',
    author_email='anascsilva@vib-ugent.be',
    url='https://www.github.com/compomics/rescore',
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
        'numpy',
        'pandas',
        'scipy',
        'scikit-learn',
        'tqdm',
        'xmltodict',
        'pyteomics'
    ],
    test_suite='ms2rescore.tests.test_ms2rescore.Tests',
    tests_require=[
        'pytest',
        'numpy',
        'pandas',
        'scipy',
        'sklearn',
        'tqdm',
        'xmltodict',
        'pyteomics',
    ],
)
