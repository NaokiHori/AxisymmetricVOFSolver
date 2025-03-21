import os
import sys

sys.path.append(os.path.abspath("./conf_params"))

project = "Axisymmetric VOF Solver"
author = "Naoki Hori"
copyright = f"2025, {author}"

from alabaster_params import html_theme
from alabaster_params import html_static_path
from alabaster_params import html_theme_options

from mathjax_params import mathjax_path
from mathjax_params import mathjax3_config

