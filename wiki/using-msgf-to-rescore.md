# Using msgf_to_rescore.py

## Necessary files

This script is written in such a way that it can call MS-GF+ by itself. To do so, the configurations file must include some additional parameters:

```
{
  "search_engine": "MSGFPlus",
  "search_engine_options": {
    "dir": "/path/to/MSGFPlus_20181015",
    "frag": "HCD",
    "path_to_modsfile": "modification.txt",
    "min_length": 8,
    "min_charge": 2,
    "max_charge": 4
  }
}
```

* `search_engine` indicates which search engine to use, and in this case should read `MSGFPlus`
* `search_engine_options` includes search settings and important paths:
   * `dir` should be the absolute path to the folder where the `MSGFPlus.jar` file is located;
   * `frag` is the fragmentation method, which can currently be either `HCD` or `CID`;
   * `path_to_modsfile` should be used if searching for modified peptides, and link an MS-GF+ modifications file;
   * `min_length`, `min_charge` and `max_charge` refers to the respective MS-GF+ search settings.

To run a search, a spectrum file in the `.MGF` format and a database in the `.fasta` format must be provided.

### I have already ran my MS-GF+ search

In this case, please open a text editor and in the file `msgf_to_rescore.py` comment the following lines like so:

```
    # logging.info("Running MS-GF+")
    # MSGF_DIR = config["search_engine_options"]["dir"]
    # run_msgfplus(MSGF_DIR, args.spec_file, args.fasta_file, config["search_engine_options"])
```

However, to execute the script, you will have to "trick" the code! Pass the following arguments as input:
* `<mgf>`: the name of your `.mzid` identifications file, but with `.mgf` as its suffix
* `<fasta>`: anything really, this file is only used when a search is ran :)

Keep in mind that you will still need an `.mgf` file to execute the ReScore pipeline!

## Execution

To generate a `.PEPREC` file that can be used to run the ReScore pipeline, execute the following command:

```
python msgf_to_rescore.py <.mgf> <.fasta> <config.json>
```

