#!/bin/bash


##
## Use a BigPurple compute node to run a Jupyter notebook and access it from your local machine.
## Can be executed through sbatch or directly.
## Run this script on the cluster to start a Jupyter notebook.
##
## Usage (direct):
## bash ./jupyter.sh
## Usage (via sbatch):
## sbatch --job-name=jupyter --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=8:00:00 ./jupyter.sh
##
## Further instructions will be printed after executing the script (check the log file if executed via sbatch).
##


#########################


# check that Jupyter notebook is available

if [ ! -x "$(command -v jupyter)" ] ; then
	echo -e "\n ERROR: 'jupyter' is not available \n" >&2
	exit 1
fi

if [ ! -x "$(command -v jupyter-notebook)" ] ; then
	echo -e "\n ERROR: 'jupyter-notebook' is not available \n" >&2
	exit 1
fi


#########################


# https://docs.ycrc.yale.edu/clusters-at-yale/guides/jupyter/
# originally modified for use on BigPurple by Paul Glick

# get tunneling info
XDG_RUNTIME_DIR=""
user=$(whoami)
node=$(hostname -s)
# port=$(shuf -i 8000-9999 -n 1)
# generate a unqiue port for each user
port=$(shuf -i 8000-9999 -n 1 --random-source <(echo "$user"))

echo -e "

Two additional steps should be perfomed on a local machine.

(1) Create an SSH tunnel in a new terminal on a local maching (there is no output):
      ssh -N -L ${port}:${node}:${port} ${user}@bigpurple.nyumc.org

(2) Access Jupyter through a web browser at:
      http://127.0.0.1:${port} (complete URL with token string will be shown below)

"

# clean up the environment and load modules or conda environments (should be a parameter)
# module purge
# module add default-environment

# classic Jupyter Notebook
# jupyter-notebook --no-browser --port=${port} --ip=${node}

# JupyterLab
jupyter lab --no-browser --port=${port} --ip=${node}


#########################



# end
