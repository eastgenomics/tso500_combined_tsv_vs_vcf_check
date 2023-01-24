# TSO500 CombinedVariantOutput tsv vs vcf check

Script to compare VCFs against TSO500 CombinedVariantOutput tsvs files. This is to find any differences in variant calls between the tsv and vcf, and output 2 separate files detailing any variants present only in either the vcf or tsv for a given run of samples.

This expects to take a path to a directory of directories as input, which each sub directory being for a run and containing all tsvs and vcfs to compare.

```
# example of dir containing 2 dirs of run files
$ ls  /tmp/tso500_reports_v1.1.0_compare/
002_230109_A01303_0138_AHKW3GDRX2_TSO500
002_230111_A01303_0139_AHKW2YDRX2_TSO500

# run comparison
$ python3 compare_tso500_tsv_to_vcf.py /tmp/tso500_reports_v1.1.0_compare/

# files output from comparing 2 runs
002_230109_A01303_0138_AHKW3GDRX2_TSO500_all_tsv_only.tsv
002_230109_A01303_0138_AHKW3GDRX2_TSO500_all_vcf_only.tsv
002_230111_A01303_0139_AHKW2YDRX2_TSO500_all_tsv_only.tsv
002_230111_A01303_0139_AHKW2YDRX2_TSO500_all_vcf_only.tsv
```