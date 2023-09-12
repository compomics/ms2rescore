from pathlib import Path

from ms2rescore.config_parser import _parse_output_path


def test__parse_output_path():
    # Ensure that test dir exists
    Path("examples/id").mkdir(parents=True, exist_ok=True)
    test_psm_file = "some/dir/psm_file.mzid"

    test_cases = [
        ("examples/id", "examples/id/psm_file.ms2rescore"),  # Existing dir
        ("examples/id/custom_stem", "examples/id/custom_stem"),  # Parent is existing dir
        ("some/other_dir", "some/other_dir/psm_file.ms2rescore"),  # None-existing dir
        (
            "some/other_dir/",
            "some/other_dir/psm_file.ms2rescore",
        ),  # None-existing dir, with trailing slash
        (None, "some/dir/psm_file.ms2rescore"),
    ]

    for output_path, expected in test_cases:
        assert _parse_output_path(output_path, test_psm_file) == expected
