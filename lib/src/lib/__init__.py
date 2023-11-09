from pathlib import Path

import yaml

path_colors = Path(__file__).parent / "colors.yaml"
with open(path_colors, "r") as f:
    colors = yaml.safe_load(f)
