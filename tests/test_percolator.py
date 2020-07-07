"""Test ms2rescore.percolator."""

from ms2rescore.percolator import PercolatorIn

TEST_CASES = [
    {
        "modified_sequences": [
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
            (
                "R.Q[-17.02655]DYC[57.02147]PMDPLGADHDDAR.L",
                "1|Gln->pyro-Glu|4|Carbamidomethyl",
                "QDYCPMDPLGADHDDAR",
            )
        ],
        "mapping": {
            42.010565: "Acetyl",
            -17.026549: "Gln->pyro-Glu",
            -18.010565: "Glu->pyro-Glu",
            57.021464: "Carbamidomethyl",
            15.994915: "Oxidation",
        }
    },
    {
        "modified_sequences": [
            ("R.LAVPIV.-", "-", "LAVPIV"),
            ("R.E[UNIMOD:27]AVPIV.-", "1|Glu->pyro-Glu", "EAVPIV"),
            ("R.IAADM[UNIMOD:35]MAISR.D", "5|Oxidation", "IAADMMAISR"),
            ("M.M[UNIMOD:1]AYGSGNK.A", "1|Acetyl", "MAYGSGNK"),
            ("-.M[UNIMOD:1]TSSEFKK.A", "1|Acetyl", "MTSSEFKK"),
            (
                "K.Q[UNIMOD:28]M[UNIMOD:35]DLC[UNIMOD:4]METIESETR.R",
                "1|Gln->pyro-Glu|2|Oxidation|5|Carbamidomethyl",
                "QMDLCMETIESETR",
            ),
            (
                "K.E[UNIMOD:27]AGDDVSM[UNIMOD:35]IARR.V",
                "1|Glu->pyro-Glu|8|Oxidation",
                "EAGDDVSMIARR",
            ),
            (
                "K.Q[UNIMOD:28]M[UNIMOD:35]FNTFNYGAANIEDGK.T",
                "1|Gln->pyro-Glu|2|Oxidation",
                "QMFNTFNYGAANIEDGK",
            ),
            (
                "-.C[UNIMOD:385]LNLQLC[UNIMOD:4]DIK.K",
                "1|Ammonia-loss|7|Carbamidomethyl",
                "CLNLQLCDIK",
            ),
            (
                "R.Q[UNIMOD:28]DYC[UNIMOD:4]PMDPLGADHDDAR.L",
                "1|Gln->pyro-Glu|4|Carbamidomethyl",
                "QDYCPMDPLGADHDDAR",
            )
        ],
        "mapping": {
            "UNIMOD:1": "Acetyl",
            "UNIMOD:28": "Gln->pyro-Glu",
            "UNIMOD:27": "Glu->pyro-Glu",
            "UNIMOD:385": "Ammonia-loss",
            "UNIMOD:4": "Carbamidomethyl",
            "UNIMOD:35": "Oxidation",
        }
    }
]


class TestPercolatorIn:
    def test_get_peprec_modifications(self):
        for case in TEST_CASES:
            pin = PercolatorIn(modification_mapping=case["mapping"])
            for mod_seq, expected_mods, _ in case["modified_sequences"]:
                assert expected_mods == pin._get_peprec_modifications(mod_seq)

    def test_get_unmodified_sequence(self):
        for case in TEST_CASES:
            pin = PercolatorIn(modification_mapping=case["mapping"])
            for mod_seq, _, expected_seq in case["modified_sequences"]:
                assert expected_seq == pin._get_unmodified_sequence(mod_seq)
