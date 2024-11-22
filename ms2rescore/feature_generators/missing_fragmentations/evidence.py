from psm_utils import Peptidoform
import numpy as np

class PeptideEvidence:
    def __init__(
            self,
            peptidoform,
            evidence_labels,
            ion_matrix=None,
            evidence=None
        ):

        if isinstance(peptidoform, str):
            self.peptidoform = Peptidoform(peptidoform)
        else:
            self.peptidoform = peptidoform

        self.evidence_labels = evidence_labels

        if evidence is not None:
            self.evidence = evidence
        else:
            self.evidence = ion_matrix.any(axis=0) 

    def __repr__(self):
        peptide_repr = ""
        amb_tag_idx = self.get_ambiguous_tag_idx()
        first_tags = [i[0] for i in amb_tag_idx]
        end_tags = [i[1]-1 for i in amb_tag_idx]

        for i, aa in enumerate(self.peptidoform.parsed_sequence):
            mod = ""
            if aa[1]:
                mod = "[" + str(aa[1][0]) + "]"
            aa_repr = aa[0]+mod

            if i in first_tags:
                peptide_repr += "<"
            peptide_repr += aa_repr
            if i in end_tags:
                peptide_repr += ">"
        peptide_repr += "/"+str(self.peptidoform.precursor_charge)
        return peptide_repr
    
    @classmethod
    def load(cls, peptide_evidence):
        
        return PeptideEvidence(
            peptidoform=peptide_evidence.peptidoform,
            evidence=peptide_evidence.evidence,
            evidence_labels=peptide_evidence.evidence_labels
        )

    def get_ambiguous_tag_idx(self, add_nterm_index=False):
        """
        Creates a list of tuples (start_index, end_index) of amino acids which can be substituted with anything with equal mass
        """
        prev_evidence = True
        isobaric_parts = []
        isobaric_part = []
        for i, aa in enumerate(self.peptidoform.parsed_sequence):

            if i == len(self.peptidoform.parsed_sequence)-1:
                break
            
            evidence = self.evidence[i]
            
            if prev_evidence:
                if not evidence:
                    isobaric_part.append(i)
            if not prev_evidence:
                if evidence:
                    isobaric_part.append(i+1)
                    isobaric_parts.append(isobaric_part)
                    isobaric_part = []
            prev_evidence=evidence

        if not prev_evidence:
            isobaric_part.append(i+1)
            isobaric_parts.append(isobaric_part)

        if add_nterm_index:
            return np.array(isobaric_parts)+1
        else:
            return np.array(isobaric_parts)
    
    @property
    def ambiguous_tags(cls):
        tags = []
        parsed_seq = cls.peptidoform.parsed_sequence

        for tag_idx in cls.get_ambiguous_tag_idx():

            tag = parsed_seq[tag_idx[0]:tag_idx[1]]

            tag_string = ""
            if tag_idx[0]==0:
                n_mod = stringify_mods(cls.peptidoform.properties["n_term"])
                if n_mod is not None:
                    tag_string += n_mod + '-'

            for aa, mod in tag:
                tag_string += f'{aa}{stringify_mods(mod)}'

            tags.append(tag_string)
        return tags

def stringify_mods(modification_list):

    s = ''
    if modification_list is None:
        return s

    for mod in modification_list:
        mod_string = f"[{mod.prefix_name}:{mod.id}]"
        s += mod_string
    return s