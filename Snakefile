#python3

## ---------------- Import Functions --------------- ## 
import sys
sys.path.append('src/')




## ---------------- Variable Setup ----------------- ##
# NOTE : these variables should be changed for any example dataset
input_vcf = 'data/test.vcf.gz'
poplabels = 'params/poplists/indiv2pop.txt'
poppanel = 'params/poplists/pop_order.txt'


## ------ Calculating Allele Frequency Table ------- ##
rule calc_af_table:
  input:
    vcf=input_vcf,
    pops = poplabels
  output:
    tmp_vcf = temp('data/test.biallelic_snps.frq.vcf.gz'),
    plink_test = temp('data/test.strat.frq.strat.gz'),
  run:
    prefix = 'data/test'
    shell('bcftools annotate --set-id \'%POS\' {input.vcf} | bcftools view -v snps -m2 -M2 | bgzip -@10 > {output.tmp_vcf}')
    shell('plink2 --vcf {output.tmp_vcf} --double-id --freq gz --within {input.pops} --out {prefix}')
    # Compute the allele frequency tables using built in functions
    

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

