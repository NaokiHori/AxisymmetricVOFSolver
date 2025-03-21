def add(macros):
    # a symbol used to denote general coordinate
    macros["gcs"] = ["{\\xi^{#1}}", 1]
    # fields
    macros["density"] = "{\\rho}"
    macros["viscosity"] = "{\\mu}"
    macros["vel"] = ["{u_{#1}}", 1]
    # scale factors
    macros["sfact"] = ["{h_{\\gcs{#1}}}", 1]
    # velocity-gradient tensor and shear-stress tensor
    macros["vgt"] = ["{l_{#1 #2}}", 2]
    macros["sst"] = ["{\\tau_{#1 #2}}", 2]
    # derivatives
    macros["pder"] = ["{\\frac{\\partial #1}{\\partial #2}}", 2]
