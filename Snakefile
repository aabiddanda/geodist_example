#python3

## ---------------- Import Functions --------------- ## 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')


# importing custom libraries and functions
import sys
sys.path.append('src/')
from reformat_plink_freq import  compute_freq_table, write_output_str
from geodist_naive import gen_geodist
from plot_geodist import GeoDistPlot


## ---------------- Variable Setup ----------------- ##
# NOTE : these variables should be changed for any example dataset
input_vcf = 'data/vcf/test.vcf.gz'
pop_labels = 'params/poplists/indiv2pop.txt'
pop_panel = 'params/poplists/pop_order.txt'

bins = "[0.0,0.05]"


## ------ Calculating Allele Frequency Table ------- ##

# TODO : print valuable information as we are going through this...
rule calc_af_table:
  input:
    vcf=input_vcf,
    pops = pop_labels,
    poporder = pop_panel
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
    gzip -d -c {input} | awk \'NR > 1 {{counts[$5]++}} END{{for(i in counts){{print i,counts[i];}}}}\' | sort -n -k1,1 | gzip > {output.geodist_cnt}
    """

## -------- Generating GeoDist Plots --------------- ##
rule plot_geodist:
  input:
    geodist_cnt = rules.count_geodist.output.geodist_cnt,
    pop_panel = pop_panel
  output:
    geodist_plot = 'plots/test.geodist.pdf'
  run:
    # NOTE : this will be changed based on the binning used...
    cur_geodist = GeoDistPlot()
    cur_geodist._add_text_data(input.geodist_cnt)
    cur_geodist._add_poplabels(input.pop_panel)
    cur_geodist._filter_data()
    cur_geodist._add_cmap()
    f, ax = plt.subplots(1,1,figsize=(4,8))
    _, nsnps, _ = cur_geodist.plot_geodist(ax)
    ax.set_yticks([0, 0.25,0.5,0.75,1.0])
    ax.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0'], fontsize=cur_geodist.y_lbl_fontsize)
    ax.set_xticklabels(cur_geodist.poplist, fontsize=cur_geodist.x_lbl_fontsize, rotation=90)
    plt.savefig(output.geodist_plot, bbox_inches='tight')

rule gen_final_plots:
  """ Rule to generate final plots for a geodist dataset """
  input:
    'plots/test.geodist.pdf'