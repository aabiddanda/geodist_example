#python3

## ---------------- Import Functions --------------- ## 
import sys
sys.path.append('src/')
from reformat_plink_freq import  compute_freq_table, write_output_str
from geodist_naive import gen_geodist


## ---------------- Variable Setup ----------------- ##
# NOTE : these variables should be changed for any example dataset
input_vcf = 'data/vcf/test.vcf.gz'
poplabels = 'params/poplists/indiv2pop.txt'
poppanel = 'params/poplists/pop_order.txt'

bins = "[0.0,0.05]"


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
rule calc_geodist_naive:
  input:
    freq_table = rules.calc_af_table.output.freq_table
  output:
    geodist_table = 'data/geodist_tbl/test.geodist.txt.gz'
  run: 
    gen_geodist(infile=input.freq_table, outfile=output.geodist_table, bins=bins)

rule count_geodist:
  input:
    rules.calc_geodist_naive.output.geodist_table
  output:
    geodist_cnt = 'data/geodist_cnts/test.geodist_cnts.txt.gz'
  shell:
    """
    gzip -d -c {input} | awk \'NR > 1 {{counts[$6]++}} END{{for(i in counts){{print i,counts[i];}}}}\' | sort -n -k1,1 | gzip > {output.geodist_cnt}
    """

rule final_count_geodist:
  input:
    'data/geodist_cnts/test.geodist_cnts.txt.gz'

## -------- Generating GeoDist Plots --------------- ##
rules plot_geodist:
  input:
    rules.count_geodist.output.geodist_cnt
  output:
    geodist_plot = ''
