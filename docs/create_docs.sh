#Automatically generate .rst files for thermal_history package, using _templates, then build the html documentation

#Required for building these docs: sphinx extensions not built in(see conf.py): sphinx read the docs theme, myst (for parsing the markdown files).
#see https://myst-parser.readthedocs.io/en/latest/sphinx/intro.html
#conda install -c conda-forge myst-parser

#and https://sphinx-rtd-theme.readthedocs.io/en/stable/installing.html
#conda install sphinx_rtd_theme


export SPHINX_APIDOC_OPTIONS="members,show-inheritance"
sphinx-apidoc -ef -t ./_templates --tocfile API_index -o ./_source ../thermal_history
make clean html
