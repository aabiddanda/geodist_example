#python3

## ---------------- Import Functions --------------- ## 
import sys
sys.path.append('src/')
from reformat_plink_freq import  compute_freq_table, write_output_str



## ---------------- Variable Setup ----------------- ##
# NOTE : these variables should be changed for any example dataset
input_vcf = 'data/vcf/test.vcf.gz'
poplabels = 'params/poplists/indiv2pop.txt'
poppanel = 'params/poplists/pop_order.txt'


## ------ Calculating Allele Frequency Table ------- ##

# TODO : print valuable information as we are going through this...
rule calc_af_table:
  input:
    vcf=input_vcf,
    pops = poplabels,
    poporder = poppanel
  output:
    tmp_vcf = temp('data/vcf/test.biallelic_snps.frq.vcf.gz'),
    plink_test = temp('data/vcf/test.frq.strat.gz'),
    freq_table = 'data/freq_table/test.freq.txt.gz'
  run:
    prefix = ''.join(input.vcf.split('.')[:-2])
    shell('bcftools annotate --set-id \'%POS\' {input.vcf} | bcftools view -v snps -m2 -M2 | bgzip -@10 > {output.tmp_vcf}')
    shell('plink2 --vcf {output.tmp_vcf} --double-id --freq gz --within {input.pops} --out {prefix}')
    # Compute the allele frequency tables using built in functions
    str_acc = compute_freq_table(freq_file = output.plink_test, pop_labels=input.poporder, sep=None, header=True)
    write_output_str(str_acc, output.freq_table)


## --- Generating Geographic Distribution Codes ---- ##
# rule calc_geodist_naive:
#   input:
#     rules.calculate_allele_frequency.output.freq_table
#   output:
#     geodist_tables

# rules count_geodist:
#     input:
#         rules.calc_geodist_naive.output.geodist_table

## -------- Generating GeoDist Plots --------------- ##
# rules plot_geodist:
#     input:
#         rules.count_geodist.output.geodist_cnt

