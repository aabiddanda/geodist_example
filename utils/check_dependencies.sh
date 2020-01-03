#!bash

#0. Checking for pip 
pip=`which pip3`

#1. Check for the list of python packages
echo "1. Checking python libraries ..."
mpl_indicator=$($pip list | grep "matplotlib")
[[ -z "$mpl_indicator" ]] && { echo "matplotlib not installed! Try \`pip install matplotlib\`"; exit 1; }
echo "matplotlib installed via pip3"

numpy_indicator=$($pip list | grep "numpy")
[[ -z "$numpy_indicator" ]] && { echo "numpy not installed! Try \`pip install numpy\`"; exit 1; }
echo "numpy installed via pip3"

snakemk_indicator=$($pip list | grep "snakemake")
[[ -z "$snakemk_indicator" ]] && { echo "snakemake not installed! Try \`pip install snakemake\`"; exit 1; }
echo "snakemake installed via pip3"

#2. Check for installed programs
echo "2. Checking installed command-line programs ..."
bcftools_indicator=$(which bcftools)
[[ -z "$bcftools_indicator" ]] && { echo "bcftools not installed! Please lookup how to install for your system"; exit 1; }
echo "bcftools ... check!"

plink_indicator=$(which plink; which plink2)
[[ -z "$plink_indicator" ]] && { echo "plink not installed! Please look up how to install plink >v1.9 for your system"; exit 1; }
echo "plink ... check!"

gzip_indicator=$(which gzip)
[[ -z "$gzip_indicator" ]] && { echo "gzip not installed! Please look up how to install plink >v1.9 for your system"; exit 1; }
echo "gzip ... check!"


awk_indicator=$(which awk)
[[ -z "$gzip_indicator" ]] && { echo "gzip not installed! Please look up how to install plink >v1.9 for your system"; exit 1; }
echo "awk ... check!" 

#3. Final message check!
echo ""
echo "All dependencies checked out. All systems go!" 
