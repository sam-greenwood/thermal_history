#!/bin/bash



echo "READING environment.yml"


env=$(grep -w 'name' environment.yml | sed -n -e 's/^.*name: //p')    #Get environment name from file

conda env create -f environment.yml                                   #Create new environment from environment.yml

source activate $env #activate


echo "INSTALLING thermal_history"

python setup.py install

#Use ipykernel to install that environment to be selectable

echo "--------------------"

echo "If you plan on using an interactive python console (e.g. jupyter notebook) you may need to install the new environment's kernel specifically. This can be done with the ipykernel python package (installed with the environment.yml file) using the command (assuming you have ipykernel installed):"

echo " "

echo "python -m ipykernel install --user --name $env --display-name \"Python ($env)\""

echo "--------------------"

echo "You can now activate the environment with:"

echo " "

echo "source activate $env"

echo "--------------------"
#kernels are saved to ~/.local/share/jupyter/kernels

#to use the correct kernel in spyder, open tools -> preferences -> python interpreter -> select path from the anaconda folder /envs/env_name/bin/python
