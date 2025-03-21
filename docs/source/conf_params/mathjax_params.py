from mathjax_macro.general import add as add_general
from mathjax_macro.momentum import add as add_momentum

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

macros = dict()

add_general(macros)
add_momentum(macros)

mathjax3_config = {"TeX": {"Macros": macros}}
