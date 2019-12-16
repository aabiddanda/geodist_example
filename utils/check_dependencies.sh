#!bash

#0. Checking for pip 
pip=`which pip`

#1. Check for the list of python packages
echo "Checking python libraries ..."
mpl_indicator=$(pip3 list | grep "matplotlib")
[[ -z "$mpl_indicator" ]] && { echo "matplotlib not installed! Try \`pip install matplotlib\`"; exit 1; }

numpy_indicator=$(pip3 list | grep "numpy")
[[ -z "$numpy_indicator" ]] && { echo "numpy not installed! Try \`pip install numpy\`"; exit 1; }

snakemk_indicator=$(pip3 list | grep "snakemake")
[[ -z "$snakemk_indicator" ]] && { echo "snakemake not installed! Try \`pip install snakemake\`"; exit 1; }

#2. Check for installed programs
echo "Checking installed programs ..."
bcftools_indicator=$(which bcftools)
[[ -z "$bcftools_indicator" ]] && { echo "bcftools not installed! Please lookup how to install for your system"; exit 1; }

plink_indicator=$(which plink; which plink2)
[[ -z "$plink_indicator" ]] && { echo "plink not installed! Try \`pip install snakemake\`"; exit 1; }


#3. Final message check!
echo "All dependencies checked out. All systems go!" 
