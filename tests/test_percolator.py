"""Test ms2rescore.percolator."""

from ms2rescore.percolator import PercolatorIn


TEST_MODIFIED_SEQUENCES = [
    ("R.LAVPIV.-", "-", "LAVPIV"),
    ("R.E[-18.01056]AVPIV.-", "1|Glu->pyro-Glu", "EAVPIV"),
    ("R.IAADM[15.99492]MAISR.D", "5|Oxidation", "IAADMMAISR"),
    ("M.M[42.01057]AYGSGNK.A", "1|Acetyl", "MAYGSGNK"),
    ("-.M[42.01057]TSSEFKK.A", "1|Acetyl", "MTSSEFKK"),
    (
        "K.Q[-17.02655]M[15.99492]DLC[57.02147]METIESETR.R",
        "1|Gln->pyro-Glu|2|Oxidation|5|Carbamidomethyl",
        "QMDLCMETIESETR",
    ),
    (
        "K.E[-18.01056]AGDDVSM[15.99492]IARR.V",
        "1|Glu->pyro-Glu|8|Oxidation",
        "EAGDDVSMIARR",
    ),
    (
        "K.Q[-17.02655]M[15.99492]FNTFNYGAANIEDGK.T",
        "1|Gln->pyro-Glu|2|Oxidation",
        "QMFNTFNYGAANIEDGK",
    ),
    (
        "-.C[-17.02655][57.02147]LNLQLC[57.02147]DIK.K",
        "0|Gln->pyro-Glu|1|Carbamidomethyl|7|Carbamidomethyl",
        "CLNLQLCDIK",
    ),
]
MODIFICATION_MAPPING = {
    42.010565: "Acetyl",
    -17.026549: "Gln->pyro-Glu",
    -18.010565: "Glu->pyro-Glu",
    57.021464: "Carbamidomethyl",
    15.994915: "Oxidation",
}


class TestPercolatorIn:
    def test_get_peprec_modifications(self):
        pin = PercolatorIn(modification_mapping=MODIFICATION_MAPPING)
        for mod_seq, expected_mods, _ in TEST_MODIFIED_SEQUENCES:
            assert expected_mods == pin.get_peprec_modifications(mod_seq)

    def test_get_unmodified_sequence(self):
        pin = PercolatorIn(modification_mapping=MODIFICATION_MAPPING)
        for mod_seq, _, expected_seq in TEST_MODIFIED_SEQUENCES:
            assert expected_seq == pin.get_unmodified_sequence(mod_seq)
