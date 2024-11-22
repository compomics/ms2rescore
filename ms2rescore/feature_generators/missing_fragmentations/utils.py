from pyteomics import mgf
import polars as pl
import os
from glob import glob
from tqdm import tqdm
import numpy as np
import multiprocessing

from rustyms import RawSpectrum, LinearPeptide, FragmentationModel

def fragments_to_polars(fragment_list, ion_types, neutral_losses, mz_array=None, intensity_array=None):
    # Can be done way quicker!

    ion_type = []
    ion_charge = []
    ion_nl = []
    ion_mass = []

    for fragment in fragment_list:
        # ignore radicals. NOTE: removing this will break annot_peaks_to_fragments
        if "·" in fragment.ion:
            continue
        elif fragment.ion == "precursor":
            # A precursor ion has only 1 ion index
            ion_type.append("p1")
        else:
            ion_type.append(fragment.ion)
        
        if fragment.neutral_loss is None:
            ion_nl.append("")
        else:
            ion_nl.append(fragment.neutral_loss)

        if mz_array is None:
            ion_mass.append(fragment.formula.mass())

        ion_charge.append(fragment.charge)

    if (mz_array is not None) and (intensity_array is not None):
        spec_polars = pl.DataFrame(
            data={
                "ion_type": [t[0]+str(c) for t, c in zip(ion_type, ion_charge)],
                "ion_index": [int(x[1:]) for x in ion_type],
                "charge": ion_charge,
                "neutral_loss": ion_nl,
                "mz": mz_array,
                "intensity": intensity_array
            },
            schema={
                "ion_type": pl.String,
                "ion_index": pl.Int8,
                "charge": pl.Int8,
                "neutral_loss": pl.String,
                "mz": pl.Float64,
                "intensity": pl.Float64
            }
        )
    else: 
        spec_polars = pl.DataFrame(
            data={
                "ion_type": [t[0]+str(c) for t, c in zip(ion_type, ion_charge)],
                "ion_index": [int(x[1:]) for x in ion_type],
                "charge": ion_charge,
                "neutral_loss": ion_nl,
                "mz": np.array(ion_mass) / np.array(ion_charge)
            },
            schema={
                "ion_type": pl.String,
                "ion_index": pl.Int8,
                "charge": pl.Int8,
                "neutral_loss": pl.String,
                "mz": pl.Float64
            }
        )
    spec_polars = spec_polars.filter(
        pl.col("ion_type").is_in(ion_types)
    )
    spec_polars = spec_polars.filter(
        pl.col("neutral_loss").is_in(neutral_losses)
    )
    return spec_polars



def get_annotated_spectrum(
        psm,
        mgf_spectrum,
        frag_model=FragmentationModel.All
    ):
    """
    Uses rustyms python bindings for peak mapping
    
    From an MGF-file containing raw MS2 spectra and a PSM, return annotated and theoretical spectra.

    Parameters
    ----------
    psm: psm_utils.PSM
        A PSM containing matching information
    mgf_spectrum: MGFIndexed
        A pyteomics mgf-formatted spectrum
    Returns
    -------
    tuple[Spectrum, Spectrum]
        - Annotated experimental spectrum
        - Theoretical fragments list
    """
    spectrum = RawSpectrum(
        title=mgf_spectrum["params"]["title"],
        num_scans=int(psm["spectrum_id"].split("=")[-1]),
        precursor_mass=mgf_spectrum["params"]["pepmass"][0],
        precursor_charge=int(mgf_spectrum["params"]["charge"][0]),
        mz_array=mgf_spectrum["m/z array"],
        intensity_array=mgf_spectrum["intensity array"],
        rt=mgf_spectrum["params"]["rtinseconds"]
    )
    peptide = LinearPeptide(psm["peptidoform"].proforma)
    return spectrum.annotate(peptide, frag_model), peptide.generate_theoretical_fragments(2, frag_model)

def annot_peaks_to_fragments(annotated_peaks):
    fragment_list = []
    mz_list = []
    intensity_list = []

    for peak in annotated_peaks:
        if len(peak.annotation) > 0:
            mz = peak.experimental_mz
            intensity = peak.intensity
            for fragment in peak.annotation:

                # ignore radicals. NOTE: removing this will break fragment_to_polars
                if "·" in fragment.ion:
                    continue

                fragment_list.append(fragment)
                mz_list.append(mz)
                intensity_list.append(intensity)
    return fragment_list, np.array(mz_list), np.array(intensity_list)

def build_array(indices, values, ion_type, n):
    # Apply a function that fills missing indices with 0.0 and constructs the array
    result = np.zeros(n)
    
    indices = np.array(indices)

    if ion_type[0] in ["x","y","z"]:
        indices = indices*(-1)
    else:
        indices = indices-1

    result[indices] = np.array(values)
    return result

def parse_to_iontype_dict(pl_df, len_pep, ion_types, value="mz"):
    """
    Parse a spectrum in polars format to a dictionary.

    Parameters
    ----------
    pl_df: polars.DataFrame
        Spectrum represented in polars dataframe. Has following columns:
        - ion_type: pl.String,
        - ion_index: pl.Int8,
        - charge: pl.Int8,
        - neutral_loss: pl.String,
        - mz: pl.Float64
        - intensity: pl.Float64
    len_pep:
        Length of the peptide annotated in the spectrum
    value:
        The column name to use for storing (either mz or intensity)
    
    Return
    ------
    dict
        - Outer keys: neutral loss
        - Inner keys: ion type
        - values: numpy arrays containing the elements sorted (0 if missing)
    """
    # Step 1: Group by 'a' and 'b' and collect the 'c' and 'd' columns as lists
    grouped_df = pl_df \
        .group_by(["neutral_loss", "ion_type"], maintain_order=True) \
        .agg([
            pl.col("ion_index").alias("indices"),
            pl.col(value).alias("values")
        ])

    # Use the 'apply' method to create the arrays and add them as a new column
    grouped_df = grouped_df \
        .with_columns([
            pl.struct(["indices", "values", "ion_type"]) \
                .map_elements(lambda x: build_array(x["indices"], x["values"], x["ion_type"], len_pep-1), return_dtype=pl.List(pl.Float64))
                .alias("array")
        ])
    # Step 3: Now pivot the DataFrame into the desired dictionary format
    # First, pivot into a dictionary structure based on 'a' and 'b' keys
    result_df = grouped_df.group_by("neutral_loss", maintain_order=True) \
        .agg([
            pl.col("ion_type"),
            pl.col("array")
        ])

    # Convert to a nested dictionary format where the outer keys are 'a' and inner keys are 'b'
    result_dict = {
        row["neutral_loss"]: {subkey: array for subkey, array in zip(row["ion_type"], row["array"])}
        for row in result_df.to_dicts()
    }
    for nl_dict in result_dict.values():
        for ion_type in ion_types:
            if ion_type not in nl_dict.keys():
                nl_dict[ion_type] = np.zeros(len_pep-1)

    return result_dict

def ion_dict_to_matrix(ion_dict, ion_types, n):
    matrix = []
    for ion_type in ion_types:
        if ion_type not in ion_dict.keys():
            matrix.append(np.zeros(n-1))
        else:
            matrix.append(ion_dict[ion_type])
    return np.array(matrix)

def matrix_to_ion_dict(matrix, ion_types):
    ion_dict = {}
    for i, ion_type in enumerate(ion_types):
        ion_dict[ion_type] = matrix[i]
    return ion_dict

def calculate_ppm(m1, m2):
    return ((m1-m2)/m1)*1e6

def mask_duplicates(matrix, preference_list):
    # Step 1: Flatten the matrix to a 1D array
    flat_matrix = matrix.flatten()

    # Step 2: Find unique values and their counts
    unique_values, counts = np.unique(flat_matrix, return_counts=True)

    # Step 3: Identify values that appear more than once
    duplicate_values = unique_values[counts > 1]

    # Step 4: Use np.isin to create a mask of locations of duplicate values
    duplicate_mask = np.isin(matrix, duplicate_values)

    # Step 5: Apply the preference list (row-based preference)
    # Create a mask to preserve only the first occurrence of duplicates based on the preference list
    final_mask = np.ones(matrix.shape, dtype=bool)  # Initialize to keep all values

    for value in duplicate_values:
        # Find all locations of the current duplicate value
        value_indices = np.argwhere(matrix == value)
        
        # Sort these locations based on the preference list
        sorted_by_preference = sorted(value_indices, key=lambda x: (preference_list.index(x[0]), x[1]))
        
        # Keep only the first occurrence based on the sorted order
        for idx in sorted_by_preference[1:]:
            final_mask[tuple(idx)] = False  # Mask out the non-preferred duplicates

    # Apply the mask to the matrix (masking out duplicates except preferred ones)
    result_matrix = np.where(final_mask, matrix, np.nan)  # You can use any placeholder like NaN for masked values
    
    return result_matrix, final_mask