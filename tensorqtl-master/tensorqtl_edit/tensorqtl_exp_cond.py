# python 3
import pandas as pd
import torch
import os
from core import *
import genotypeio, eigenmt
import cis_ase
import cis_exon
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='linear regression: Z-score ~ Heter + COV1+ ... + COVn\n using tensorqtl')
parser.add_argument('--mode', '-m', choices=['cis', 'cis_nominal', 'cis_independent', 'trans'],
                    help='mode: cis cis_nominal cis_independent', required=True)
parser.add_argument('--input', '-i', help='exon_id z-value bed file ', action='store')
parser.add_argument('--geno', '-g', help='genoytpe file bfile prefix')
parser.add_argument('--cov', '-c', help='cov file ')
parser.add_argument('--out', '-o', help='output perfix file ', action='store')
parser.add_argument('--output_dir', default='.', help='Output directory')
# parser.add_argument('--maf', '-m', help='MAF filter remove low maf', default='0.001')
parser.add_argument('--chr', '-r', help='select chr to run (e.g. chr1 (default All)', default='All')
parser.add_argument('--cis_output', help='cis output results with FDR (only for mode: indep', default=None)
parser.add_argument('--cond_exp', help='gene expression bed.file as condition', default=None)
parser.add_argument('--fdr', help='fdr threshold (only for mode: indep), default=0.05', default=0.05)
parser.add_argument('--p_beta_th', help='p_beta threshold genome wide @ FDR (only for mode: indep)', default=None)
parser.add_argument('--seed', default=None, type=int, help='Seed for permutations.')
args = parser.parse_args()

all_chrs_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

# logger = SimpleLogger()

expression_bed = args.input
covariates_file = args.cov
plink_prefix_path = args.geno
prefix = args.out
chr_id = args.chr
mode = args.mode
in_cis_addr = args.cis_output
fdr = float(args.fdr)

# check output_dir
if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)
logger = SimpleLogger(os.path.join(args.output_dir, prefix+'.tensorQTL.{}.log'.format(args.mode)))


if args.chr == 'All':
    excluded_chr_list = None
else:
    all_chrs_list.remove(chr_id)
    excluded_chr_list = all_chrs_list

logger.write('[{}] Running TensorQTL: {}-QTL mapping'.format(datetime.now().strftime("%b %d %H:%M:%S"), args.mode.split('_')[0]))


# load phenotypes and covariates
logger.write('  * reading phenotypes ({})'.format(args.input))
phenotype_df, phenotype_pos_df = read_phenotype_bed(args.input)
# covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

if args.cov is not None:
    logger.write('  * reading covariates ({})'.format(args.cov))
    covariates_df = pd.read_csv(args.cov, sep='\t', index_col=0).T
    assert np.all(phenotype_df.columns == covariates_df.index)


if args.cond_exp is not None:
    express_df, express_pos_df = read_phenotype_bed(args.cond_exp)
    assert np.all(phenotype_df.columns == express_df.columns)

pr = genotypeio.PlinkReader(plink_prefix_path, exclude_chrs=excluded_chr_list)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

if args.mode == 'cis_normal':
    if args.cis_output and os.path.exists(args.cis_output):  # only run on the significant phenotype
        fdr = args.fdr
        logger.write(' * only using FDR < {} phenotype'.format(fdr))
        cis_df = pd.read_csv(args.cis_output, sep='\t', index_col=0)
        fdr_col = 'qval'
        signif_df = cis_df[cis_df[fdr_col] <= fdr].copy()
        ix = phenotype_df.index[phenotype_df.index.isin(signif_df.index)]
        logger.write('  * {}/{} significant phenotypes'.format(signif_df.shape[0], cis_df.shape[0]))
        phenotype_df = phenotype_df.loc[ix]
        phenotype_pos_df = phenotype_pos_df.loc[ix]

        cis_exon.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
                            covariates_df=covariates_df,
                            group_s=None,  run_eigenmt=True, output_dir=args.output_dir,
                            write_top=True, write_stats=not args.best_only, logger=logger, verbose=True,
                            express_df=express_df, express_pos_df=express_pos_df)
else:
    print('Please use other script!')