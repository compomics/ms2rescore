#! python

from setuptools import setup

setup(
    name='ms2rescore',
    version='1.0',
    description='MS²ReScore: Search engine independent rescoring of PSMs by using predicted MS² peak intensities.',
    author='Ana Sílvia C. Silva',
    author_email='anascsilva@vib-ugent.be',
    url='https://www.github.com/compomics/rescore',
    packages=['ms2rescore'],
    scripts=['ms2rescore/ms2rescore'],
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta"
    ],
    test_suite='ms2rescore.tests.test_ms2rescore.Tests',
    tests_require=['pytest'],
)
