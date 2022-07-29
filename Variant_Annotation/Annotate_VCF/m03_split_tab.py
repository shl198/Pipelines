import gzip
import argparse,glob
import pandas as pd 


parser = argparse.ArgumentParser(description='split code, regulatory and other variants')

parser.add_argument('-i','--input',action='store',dest='tab',help='input tab file')
parser.add_argument('-c','--code',action='store',dest='code_table',help='output code table file')
parser.add_argument('-r','--regu',action='store',dest='regu_table',help='output regu table file')
parser.add_argument('-n','--ncode',action='store',dest='ncode_table',help='output non-code table file')

args = parser.parse_args()

tab_file = args.tab
code_tab_file = args.code_table
regu_tab_file = args.regu_table
ncode_tab_file = args.ncode_table

coding_effect = ['transcript_ablation','splice_acceptor_variant',
    'splice_donor_variant','stop_gained',
    'splice_region_variant','frameshift_variant',
    'stop_lost','start_lost',
    'transcript_amplification',
    'inframe_insertion','inframe_deletion',
    'missense_variant','protein_altering_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant','stop_retained_variant',
    'synonymous_variant','coding_sequence_variant']

regulate_effect = ['mature_miRNA_variant',
    '5_prime_UTR_variant','3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'TFBS_ablation','TFBS_amplification',
    'TF_binding_site_variant','regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation']

inter_intron_effect = ['intron_variant','intergenic_variant',
               'downstream_gene_variant','upstream_gene_variant']

df = pd.read_csv(tab_file, sep='\t',header=0,compression='gzip')
code_df = df.query('Consequence in @coding_effect')
code_df.to_csv(code_tab_file,sep='\t',index=False,compression='gzip')
regu_df = df.query('Consequence in @regulate_effect')
regu_df.to_csv(regu_tab_file,sep='\t',index=False,compression='gzip')
ncode_df = df.query('Consequence in @inter_intron_effect')
ncode_df.to_csv(ncode_tab_file,sep='\t',index=False,compression='gzip')
