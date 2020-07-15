"""Convert PIN file to PEPREC for MS²ReScore."""

import logging
from typing import Dict, Tuple

import click

from ms2rescore.percolator import PercolatorIn


SPEC_ID_PATTERNS = {
    "generic": r".+_([0-9]+)_[0-9]+_[0-9]+",
    "tandem": r".+_([0-9]+)_[0-9]+_[0-9]+",
    "msgfplus": r".+_SII_([0-9]+)_[0-9]+_[0-9]+_[0-9]+",
}


def parse_mgf(path_to_mgf: str) -> Tuple[Dict[int, str], Dict[int, float]]:
    """Parse MGF file to extract title and retention time fields, by spectrum index."""
    titles = dict()
    retention_times = dict()
    with open(path_to_mgf, "rt") as mgf_in:
        index = 0
        for line in mgf_in:
            if line[0] == "B":
                if line.strip() == "BEGIN IONS":
                    index += 1
            if line[0] == "T":
                if line.startswith("TITLE="):
                    titles[index] = line[6:].strip()
            if line[0] == "R":
                if line.startswith("RTINSECONDS="):
                    retention_times[index] = float(line[12:].strip())
    return titles, retention_times


def pipeline(config: Dict) -> Tuple[str, str]:
    """Convert PIN to PEPREC for MS²ReScore."""
    logging.info("Parsing PIN file...")
    # PIN to PEPREC
    pin = PercolatorIn(config["general"]["identification_file"])
    peprec = pin.to_peptide_record(spectrum_index_pattern=SPEC_ID_PATTERNS["generic"])

    # Map MGF titles and observed retention times
    titles, retention_times = parse_mgf(config["general"]["mgf_file"])
    peprec.df["observed_retention_time"] = peprec.df["spec_id"].map(retention_times)
    peprec.df["spec_id"] = peprec.df["spec_id"].map(titles)
    assert (
        ~peprec.df["observed_retention_time"].isna().any()
    ), "Could not map all MGF retention times to spectrum indices."
    assert (
        ~peprec.df["spec_id"].isna().any()
    ), "Could not map all MGF titles to spectrum indices."

    # Save PEPREC
    peprec_path = config["general"]["output_filename"] + ".peprec"
    peprec.to_csv(peprec_path)

    return peprec_path, config["general"]["mgf_file"]


@click.command()
@click.argument("input-pin")
@click.argument("output-peprec")
def main(**kwargs):
    """Convert PIN to PEPREC."""
    pin = PercolatorIn(kwargs["input_pin"])
    peprec = pin.to_peptide_record()
    peprec.to_csv(kwargs["output_peprec"])


if __name__ == "__main__":
    main()
