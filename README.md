# geodist_example
-------------------------------

An example pipeline to generate "GeoDist"-style plots from a VCF

## Requirements

### Python Packages

The following packages must be installed via `pip` or `pip3` 

1. `snakemake`
2. `numpy`
3. `matplotlib`

### External Software

Parts of the GeoDist pipeline use the following tools 

1. `bcftools`
2. `plink` (> v1.9)
3. `gzip`
4. `awk`

If installing on Mac OSX and using [homebrew](https://brew.sh/) it is possible to install all of the programs using:

```
brew install bcftools plink2 gzip awk
```

On a Linux machine the programs can be installed using `apt-get`.

### Checking Dependencies

If you are interested in checking whether your system is able to run this pipeline, 
we have provided a shell script to check whether the python libraries and command line programs are installed and visible in your path.
```
bash utils/check_dependencies.sh
```

Note : this script uses `pip3` in order to check for the `python` modules

## Running the Example Pipeline

Once all of the above dependencies have been installed, run:

```
snakemake test_final --dryrun
```

This should list all of the steps in the pipeline to generate intermediate files and the eventual PDF plot. To actually run the pipeline, remove the `--dryrun` flag above. 

After the pipeline has run, you should see the file `plots/test.geodist.pdf` as a resulting output. It should finish within 30 seconds. 

## File Descriptions

### Input Files 

#### VCF 

The input VCF file should contain all of the individuals and variation that a user may want to plot. For more information on the VCF file specification please see the [official documentation]()

#### Population Labels

The example file `params/poplists/indiv2pop.txt` is a three-column file that contains a single row for each individual of the VCF file

#### Population Panel (Population Order)

The example file `params/poplists/pop_order.txt` is a one-column file that shows an ordering for each of the superpopulation labels contained in the third column of the previously described population label file. 

### Intermediate Files

#### Frequency Table

The allele frequency table (`data/freq_table/test.freq.txt.gz` after running the pipeline)

#### GeoDist Table


#### GeoDist Counts



## Contact

If you have any questions running this pipeline please submit an issue or contact Arjun Biddanda 