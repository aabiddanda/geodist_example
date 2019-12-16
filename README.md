# geodist_example
-------------------------------

An example pipeline to generate "GeoDist"-style plots from a VCF

## Requirements

### Python Packages

The following packages may be installed via `pip` 
- snakemake
- numpy
- matplotlib

### External Software

Many parts of the GeoDist pipeline use the following tools 

1. `bcftools`
2. `plink` (> v1.9)
3. `gzip`
4. `awk`

### Checking Dependencies

If you are interested in checking whether your system is able to run this pipeline, 
we have provided a shell script to check whether the python libraries and command line programs are installed and visible in a path. 

## Contact

If you have any questions running this pipeline please submit an issue or contact Arjun Biddanda 