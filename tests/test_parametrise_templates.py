import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml
from deepdiff import DeepDiff

from lib.template import parametrise_template

EXPECTED_DIR = Path(__file__).parent / "_files"
TEMP_DIR = Path(__file__).parent / "temp"
TEMPLATE_DIR = Path(__file__).parent.parent / "data" / "templates"


def compare_yaml_files(path_file1, path_file2):
    with open(path_file1, "r") as file1, open(path_file2, "r") as file2:
        content1 = yaml.safe_load(file1)
        content2 = yaml.safe_load(file2)

        diff = DeepDiff(content1, content2, ignore_order=True)

        return diff


@pytest.fixture(autouse=True)
def prepare_temp_dir():
    """Clean up TEMP_DIR"""
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)

    TEMP_DIR.mkdir()


def test_parametrise_templates():
    templates = [
        "locations-wind-onshore.yaml",
        "locations-wind-offshore-deep.yaml",
        "locations-wind-offshore-shallow.yaml",
    ]

    files_potentials = [
        "potentials-wind-onshore.csv",
        "potentials-wind-offshore-deep.csv",
        "potentials-wind-offshore-shallow.csv",
    ]

    for template, file_potentials in zip(templates, files_potentials):
        potentials = pd.read_csv(EXPECTED_DIR / file_potentials, index_col=0)

        parametrise_template(
            TEMPLATE_DIR / template,
            TEMP_DIR / template,
            potentials=potentials,
        )

        produced = TEMP_DIR / template
        expected = EXPECTED_DIR / template

        diff = compare_yaml_files(produced, expected)

        assert not diff
