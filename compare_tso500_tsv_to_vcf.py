"""
Script to compare CombinedVariantOutput tsv files against compressed VCFs to
identify variants unique to either.

This takes a directory of directories as input (one per run), containing all
vcfs and tsvs together.

Will output up to 3 files:
    - all_tsv_only.tsv => original entries of variants only found in the tsvs
    - all_vcf_only.tsv => original entries of variants only found in the vcfs
    - all_tsv_no_annotation.tsv => original entires of variants in the tsvs
        for which there is no annotation and are not compared
"""
import gzip
import os
from pathlib import Path
import sys

import pandas as pd


def read_tsv(tsv):
    """
    Read in CombinedVariantOutput.tsv to a dataframe

    Parameters
    ----------
    tsv : str
        path to tsv file to read in
    
    Returns
    -------
    pd.DataFrame
        dataframe of SNVs read in from CombinedVariantOutput tsv
    """
    variants = []

    with open(tsv) as fh:
        # get just the variants from the tsv after the line '[Small Variants]'
        snvs = False
        for line in fh.readlines():
            if snvs:
                variants.append(line)
            if line.startswith('[Small Variants]'):
                snvs = True

    variants = [x.split('\t') for x in variants if set(x) != {'\t', '\n'}]
    variants = variants[1:]

    # if no variants found this will be returned as [['NA', '', '\n']] =>
    # set it to an empty list
    if variants == [['NA', '', '\n']]:
        variants = []

    columns = [
        'Gene', 'CHROM', 'POS', 'REF', 'ALT', 'Allele Frequency', 'Depth',
        'P-Dot Notation', 'C-Dot Notation', 'Consequence(s)', 'Affected Exon(s)'
    ]

    df = pd.DataFrame(variants, columns=columns)

    df['CHROM'] = df['CHROM'].apply(lambda x: x.replace('chr', ''))
    df['POS'] = pd.to_numeric(df['POS'])

    df['Affected Exon(s)'] = df['Affected Exon(s)'].str.replace('\n', '')

    return df


def read_vcf(vcf):
    """
    Read in compressed VCF to a dataframe

    Parameters
    ----------
    vcf : str
        path to vcf file to read in

    Returns
    -------
    pd.DataFrame
        dataframe of variants read from vcf
    """
    with gzip.open(vcf) as fh:
        # get column names
        for line in fh.readlines():
            line = line.decode()
            if line.startswith('#CHROM'):
                column_names = line.strip('#').split('\t')
                column_names[-1] = 'SAMPLE_FORMAT'
                break

    df = pd.read_csv(
        vcf, sep='\t', comment='#', names=column_names, compression='infer'
    )

    df['CHROM'] = df['CHROM'].astype(str)

    return df


def match_tsvs_and_vcfs(files):
    """
    Match tsvs and vcfs for samples from list of files in directory.

    Works based off prefix of filename of both tsv and vcf

    Parameters
    ----------
    files : list
        glob list of tsv and vcf files from a given dir

    Returns
    -------
    tsvs : list
        list of tsv files
    vcfs : list
        list of vcf files
    """
    tsvs = sorted([x for x in files if x.endswith('.tsv')])
    vcfs = sorted([x for x in files if x.endswith('.vcf.gz')])

    # we expect tsvs to be named as SAMPLE1_CominedVariantOutput.tsv and
    # vcfs to be named SAMPLE1-more-fields.vcf.gz, therefore match on the
    # sample ID from both
    tsv_prefixes = [Path(x).name.split('_')[0] for x in tsvs]
    vcf_prefixes = [Path(x).name.split('-')[0] for x in vcfs]


    # get the samples that have both tsvs and vcfs to compare
    common = list(set(tsv_prefixes) & set(vcf_prefixes))
    tsvs = [x for x in tsvs if Path(x).name.split('_')[0] in common]
    vcfs = [x for x in vcfs if Path(x).name.split('-')[0] in common]

    return tsvs, vcfs


def get_mismatch_variants(tsv_df, vcf_df):
    """
    Get variants that are mismatched (i.e. unique to tsv or vcf)

    Parameters
    ----------
    tsv_df : pd.DataFrame
        df of variants from tsv
    vcf_df : pd.DataFrame
        df of variants from vcf

    Returns
    -------
    tsv_only : pd.DataFrame
        df of variants present only in the tsv (i.e. missing from the vcf)
    vcf_only : pd.DataFrame
        df of variants present only in the vcf (i.e. missing from the tsv)
    """
    # get variants present only in tsv and vcf
    merge_cols = ['CHROM', 'POS', 'REF', 'ALT']
    tsv_only = tsv_df.merge(
        vcf_df[merge_cols], on=merge_cols, how='outer', indicator=True
        ).loc[lambda x: x.pop('_merge').eq('left_only')]

    vcf_only = vcf_df.merge(
        tsv_df[merge_cols], on=merge_cols, how='outer', indicator=True
        ).loc[lambda x: x.pop('_merge').eq('left_only')]

    return tsv_only, vcf_only


def main():

    # should be a directory of individual run directories
    all_runs_dir = Path(sys.argv[1]).absolute()

    # counters for printing at the end
    samples = 0
    runs = 0
    all_tsv_issues = 0
    all_vcf_issues = 0
    all_no_annotation_issues = 0

    for run_dir in os.listdir(all_runs_dir):
        run_dir = os.path.join(all_runs_dir, run_dir)
        run_name = run_dir.split('/')[-1]
        if not os.path.isdir(run_dir):
            # in case of any bonus files that aren't directories
            continue

        print(f"\nChecking run {run_name} ({run_dir})")

        files = [os.path.join(run_dir, x) for x in os.listdir(run_dir)]

        tsvs, vcfs = match_tsvs_and_vcfs(files)

        # empty dfs to add all mismatches to
        all_tsv_only = pd.DataFrame(
            columns=[
                'Sample', 'Gene', 'CHROM', 'POS', 'REF', 'ALT',
                'Allele Frequency', 'Depth', 'P-Dot Notation',
                'C-Dot Notation', 'Consequence(s)', 'Affected Exon(s)'
            ]
        )

        all_vcf_only = pd.DataFrame(
            columns=[
                'SAMPLE', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                'FILTER', 'INFO', 'FORMAT', 'SAMPLE_FORMAT'
            ]
        )

        # some variants in CombinedVariantOutput files have no real annotation
        # add these to their own df to write to a separate file
        all_tsv_no_annotation = pd.DataFrame(
            columns=[
                'Sample', 'Gene', 'CHROM', 'POS', 'REF', 'ALT',
                'Allele Frequency', 'Depth', 'P-Dot Notation',
                'C-Dot Notation', 'Consequence(s)', 'Affected Exon(s)'
            ]
        )

        for tsv, vcf in zip(tsvs, vcfs):
            # sense check tsv and vcf for same sample
            if Path(tsv).name.split('_')[0] != Path(vcf).name.split('-')[0]:
                print(
                    f"Error: tsv and vcf file prefixes do not match:\n{vcf}\n{tsv}"
                )
                continue

            tsv_df = read_tsv(tsv)
            vcf_df = read_vcf(vcf)

            tsv_only, vcf_only = get_mismatch_variants(tsv_df, vcf_df)

            # get any rows only in VCF AND not rescued in rescue app as
            # we know these will only be in the vcf since we are rescuing them
            # v1.0.0 reports workflow => tagged 'OPA'
            # v1.1.0 reports workflow => tagged 'rescued'
            vcf_only = vcf_only[~vcf_only['FILTER'].str.contains('OPA')]
            vcf_only = vcf_only[~vcf_only['FILTER'].str.contains('rescued')]

            # remove variants with no annotation => weird Illumina variants
            # add these to their own dataframe to dump out at the end
            no_annotation_only = tsv_only[tsv_only['Gene'] == '']
            tsv_only = tsv_only[tsv_only['Gene'] != '']

            if len(no_annotation_only.index) > 0:
                # some variants from tsv with no annotation, add to df with
                # sample ID as first column
                sample_col = [Path(tsv).name.split('_')[0]] * len(no_annotation_only.index)
                no_annotation_only.insert(0, 'Sample', sample_col)
                all_tsv_no_annotation = pd.concat(
                    [all_tsv_no_annotation, no_annotation_only])


            if len(tsv_only.index) > 0:
                # some variants only in the tsv, add to df with
                # sample ID as first column
                sample_col = [Path(tsv).name.split('_')[0]] * len(tsv_only.index)
                tsv_only.insert(0, 'Sample', sample_col)
                all_tsv_only = pd.concat([all_tsv_only, tsv_only])

            if len(vcf_only.index) > 0:
                # some variants only in our vcf
                sample_col = [Path(tsv).name.split('_')[0]] * len(vcf_only.index)
                vcf_only.insert(0, 'SAMPLE', sample_col)
                all_vcf_only = pd.concat([all_vcf_only, vcf_only])

            samples += 1

        runs += 1

        print(f"\nTotal tsv mismatch for run {run_name}: {len(all_tsv_only.index)}")
        print(f"Total vcf mismatch for run {run_name}: {len(all_vcf_only.index)}\n\n")

        if len(all_tsv_only.index) > 0:
            all_tsv_issues += len(all_tsv_only.index)

            all_tsv_only.sort_values(by=['Gene', 'POS'], inplace=True)

            all_tsv_only.to_csv(
                f'{run_name}_all_tsv_only.tsv', mode='w',
                sep='\t', index=False, header=False
            )

        if len(all_vcf_only.index) > 0:
            all_vcf_issues += len(all_vcf_only.index)
            all_vcf_only.to_csv(
                f'{run_name}_all_vcf_only.tsv', mode='w',
                sep='\t', index=False, header=False
            )

        if len(all_tsv_no_annotation.index) > 0:
            all_no_annotation_issues += len(all_tsv_no_annotation.index)
            all_tsv_no_annotation.to_csv(
                f'{run_name}_all_tsv_no_annotation.tsv', mode='w',
                sep='\t', index=False, header=False
            )

    print(f"Total samples checked: {samples}")
    print(f"Total runs checked: {runs}")
    print(f"Total tsv mismatch: {all_tsv_issues}")
    print(f"Total vcf mismatch: {all_vcf_issues}")
    print(f"Total tsv variants w/ no annotation: {all_no_annotation_issues}")

if __name__ == "__main__":
    main()
