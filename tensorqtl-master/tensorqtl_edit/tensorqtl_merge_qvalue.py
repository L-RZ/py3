import pandas as pd
import torch
import post

import argparse

parser = argparse.ArgumentParser(description='merge results to one result and calculate q-value (FDR)')

parser.add_argument('--input', '-i', help='cis p-value results (e.g. -i res1.txt res2.txt ...',
                    action='store', nargs='+')
parser.add_argument('--output', '-o', help='output file address')
parser.add_argument('--fdr', default=0.05, type=float,  help='FDR, default=0.05')
parser.add_argument('--qvalue_lambda',  type=float, default=None,
                    help='lambda parameter for pi0est in qvalue. defulat=None (GTEX:0.85)')

args = parser.parse_args()
in_file_l = args.input
fdr = float(args.fdr)
qvalue_lambda = args.qvalue_lambda
if qvalue_lambda is not None:
    qvalue_lambda = float(args.qvalue_lambda)
out_addr = args.output

df_list = [pd.read_csv(each_file, sep='\t', index_col=0) for each_file in in_file_l]
cis_df = pd.concat(df_list)

post.calculate_qvalues(cis_df, qvalue_lambda=qvalue_lambda, fdr=fdr)

cis_df.to_csv(out_addr, sep='\t')
