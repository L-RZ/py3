# python 3
import pandas as pd
import torch
import os
from core import *
import genotypeio, eigenmt
import cis_ase
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='linear regression: Z-score ~ Heter + COV1+ ... + COVn\n using tensorqtl')
parser.add_argument('--input', '-i', help='exon_id z-value bed file ', action='store')
parser.add_argument('--geno', '-g', help='genoytpe file bfile prefix')
parser.add_argument('--cov', '-c', help='cov file ')
parser.add_argument('--out', '-o', help='output perfix file ', action='store')
parser.add_argument('--output_dir', default='.', help='Output directory')
# parser.add_argument('--maf', '-m', help='MAF filter remove low maf', default='0.001')
parser.add_argument('--chr', '-r', help='select chr to run (e.g. chr1 (default All)', default='All')
parser.add_argument('--mode', '-m', choices=['cis', 'cis_nominal', 'cis_independent', 'trans'],
                    help='mode: cis cis_nominal cis_independent', required=True)
parser.add_argument('--cis_output', help='cis output results with FDR (only for mode: indep', default=None)
parser.add_argument('--fdr', help='fdr threshold (only for mode: indep), default=0.05', default=0.05)
parser.add_argument('--p_beta_th', help='p_beta threshold genome wide @ FDR (only for mode: indep)', default=None)
parser.add_argument('--seed', default=None, type=int, help='Seed for permutations.')
args = parser.parse_args()

all_chrs_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

logger = SimpleLogger()

expression_bed = args.input
covariates_file = args.cov
plink_prefix_path = args.geno
prefix = args.out
chr_id = args.chr
mode = args.mode
in_cis_addr = args.cis_output
fdr = float(args.fdr)

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

pr = genotypeio.PlinkReader(plink_prefix_path, exclude_chrs=excluded_chr_list)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

if mode == 'cis':
    # cis-QTL: empirical p-values for phenotypes
    if excluded_chr_list:
        cis_df = cis_ase.map_cis(genotype_df, variant_df,
                             phenotype_df.loc[phenotype_pos_df['chr'] == chr_id],
                             phenotype_pos_df.loc[phenotype_pos_df['chr'] == chr_id],
                             covariates_df=covariates_df, seed=args.seed
                             )
    else:
        cis_df = cis_ase.map_cis(genotype_df, variant_df,
                             phenotype_df,
                             phenotype_pos_df,
                             covariates_df=covariates_df, seed=args.seed
                             )
    out_file = os.path.join(args.output_dir, prefix + '.cis_qtl.txt.gz')
    cis_df.to_csv(out_file, sep='\t')
elif mode == 'cis_independent':
    cis_df = pd.read_csv(in_cis_addr, sep='\t', index_col=0)
    if excluded_chr_list:
        if args.p_beta_th:
            p_beta_th = float(args.p_beta_th)
            indep_df = cis_ase.map_independent(genotype_df, variant_df, cis_df,
                                           phenotype_df.loc[phenotype_pos_df['chr'] == chr_id],
                                           phenotype_pos_df.loc[phenotype_pos_df['chr'] == chr_id],
                                           covariates_df=covariates_df, fdr=fdr, signif_th=p_beta_th,
                                           seed=args.seed
                                           )
        else:
            logger.write('*  waring: p_beta_th from Chr / region instead of Genome-wide')
            indep_df = cis_ase.map_independent(genotype_df, variant_df, cis_df,
                                           phenotype_df.loc[phenotype_pos_df['chr'] == chr_id],
                                           phenotype_pos_df.loc[phenotype_pos_df['chr'] == chr_id],
                                           covariates_df=covariates_df, fdr=fdr, seed=args.seed
                                           )


    else:
        indep_df = cis_ase.map_independent(genotype_df, variant_df, cis_df,
                                       phenotype_df, phenotype_pos_df,
                                       covariates_df=covariates_df, fdr=fdr, seed=args.seed
                                       )
    out_file = os.path.join(args.output_dir, prefix + '.cis_independent_qtl.txt.gz')
    indep_df.to_csv(out_file, sep='\t', index=False)

logger.write('[{}] Finished mapping'.format(datetime.now().strftime("%b %d %H:%M:%S")))