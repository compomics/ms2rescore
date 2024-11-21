# Reproducing manuscript results

The manuscript titled _Accurate peptide fragmentation predictions allow data driven approaches to replace and improve upon proteomics search engine scoring functions_ shows the results obtained with this tool on two datasets: the _Pyrococcus furiosus_ standard, and a real-world dataset. In this section you'll find instructions to reproduce those results. To visualize them (or any other files processed with this pipeline), feel free to take advantage of the [Jupyter Notebook](https://github.com/anasilviacs/rescore/blob/pub/manuscript/reproducing-results.ipynb) included in this repository.

## _Pyrococcus furiosus_

The spectrum file for this dataset can be downloaded from PRIDE, under identifier [PXD001077](https://www.ebi.ac.uk/pride/archive/projects/PXD001077). Download the file `Velos005137.mgf`.

Choose any search engine of your preference to search these spectra against a database. We present results from a search done against all sequences from this organism (reviewed and unreviewed) obtained from [UniProt](https://www.uniprot.org/uniprot/?query=pyrfu&fil=organism%3A%22Pyrococcus+furiosus+%28strain+ATCC+43587+%2F+DSM+3638+%2F+JCM+8422+%2F+Vc1%29+%5B186497%5D%22&sort=score).

Our search was done with [MS-GF+](https://omics.pnl.gov/software/ms-gf), allowing for carbamidomethylation and oxidation (maximum 2). From the resulting `mzid` file, a list of PSMs can be extracted. The information needed to proceed is as follows:

- Spectrum identifier (corresponding to the `TITLE` field in the peak file)
- Peptide sequence (without trailing aminoacids)
- Peptide charge
- Modifications

The PEPREC format corresponding to our search is available in the folder `manuscript/pyrfu` from this repository, under the name `Velos005137.PEPREC`. Please refer to [ms2pip_c](https://github.com/anasilviacs/ms2pip_c/tree/6f037dc2d0797cd25061aaed8091d625123971e1) for more information of the PEPREC format.

Following this the configuration file should be adapted to the current goal. If the installation instructions were followed, the field `dir` should be correct. In this case the spectra are fragmented in HCD so the `frag` field should also be correct. Please change the `num_cpu` field to the number of cores your machine has available, in order to optimize the parallel routines existent throughout the code.

With these two files (mgf and PEPREC), the pipeline can be executed by issuing the following command:

```
python ../../driver.py Velos005137.mgf Velos005137.PEPREC config.json
```

This results in a file named `Velos005137.pin`, which is the input for Percolator. This file contains peptide and spectrum information and identifiers, along with the spectral similarity features calculated. To use Percolator, issue the following command:

```
percolator Velos005137.pin --trainFDR 0.01 -m Velos005137.pout -M Velos005137.pout_dec -w Velos005137.weights
```

The resulting files contain target and decoy PSMs, as well as the weights the model attributed to each feature in the feature matrix.

## HEK sample

The spectra corresponding to this experiment can be found in PRIDE, under identifier [PXD001468](https://www.ebi.ac.uk/pride/archive/projects/PXD001468). Here the spectra are split into several mgf files; the results shown in the manuscript were obtained by processing them all together. This task is significantly more computationally intensive so it does take some time. If this isn't and issue, we make the concatenated spectrum file available from our servers: [genesis.ugent.be/uvpublicdata/silvia/PXD001468/all.mgf](http://genesis.ugent.be/uvpublicdata/silvia/PXD001468/all.mgf). In case this task proves too computationally demanding, please process each (or only one) spectrum file individually.

The file `all.PEPREC` (which can be downloaded from [genesis.ugent.be/uvpublicdata/silvia/PXD001468/all.PEPREC](http://genesis.ugent.be/uvpublicdata/silvia/PXD001468/all.PEPREC)), was compiled from a search done with [MS-GF+](https://omics.pnl.gov/software/ms-gf), allowing for carbamidomethylation and oxidation (maximum 2). The configuration file in this folder should not need to be adjusted, with the exception of the `num_cpu` parameter which should reflect the number of CPU's you can dedicate to this task. Keep in mind that for the file including all the spectra there are > 1000000 PSMs and spectrum predictions and similarity metric computations will take a long time.

Issue the following command:

```
python ../../driver.py all.mgf all.PEPREC config.json
```

The result should be a file named `all.pin`. To run Percolator in an efficient way, you can tune its parameters. We ran it as follows:

```
percolator all.pin --trainFDR 0.01 --subset-max-train 200000 -m all.pout -M all.pout_dec -w all.weights
```

The resulting files contain target and decoy PSMs, as well as the weights the model attributed to each feature in the feature matrix.
