# TSO500 CombinedVariantOutput tsv vs vcf check

Script to compare VCFs against TSO500 CombinedVariantOutput tsvs files. This is to find any differences in variant calls between the tsv and vcf, and output 2 separate files detailing any variants present only in either the vcf or tsv for a given run of samples.

This expects to take a path to a directory of directories as input, which each sub directory being for a run and containing all tsvs and vcfs to compare.

```
# example of dir containing dirs of run files
$ tree /tmp/tso500_reports_v1.1.0_compare/
/tmp/tso500_reports_v1.1.0_compare/
└── 002_230109_A01303_0138_AHKW3GDRX2_TSO500
    ├── HD753-216776-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── HD753_CombinedVariantOutput.tsv
    ├── MB5081E-216768-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5081E_CombinedVariantOutput.tsv
    ├── MB5255-216769-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5255_CombinedVariantOutput.tsv
    ├── MB5269-216770-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5269_CombinedVariantOutput.tsv
    ├── MB5270-216771-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5270_CombinedVariantOutput.tsv
    ├── MB5273-216772-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5273_CombinedVariantOutput.tsv
    ├── MB5275-216773-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5275_CombinedVariantOutput.tsv
    ├── MB5285-216774-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    ├── MB5285_CombinedVariantOutput.tsv
    ├── MB5286-216775-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.gz
    └── MB5286_CombinedVariantOutput.tsv
└── 002_230119_A01303_0142_BHVLMWDRX2_TSO500
    ├── ....


# run comparison script
$ python3 compare_tso500_tsv_to_vcf.py /tmp/tso500_reports_v1.1.0_compare/

# files output from comparing 2 runs
{run_id}_all_tsv_only.tsv
{run_id}_all_vcf_only.tsv
{run_id}_all_tsv_no_annotation.tsv
```

### Output files
The script outputs 3 files per run checked:
* `_all_tsv_only.tsv` - all variants present only in the CombinedVariantOutput.tsv files (i.e. missing from the vcf), this will have one line per sample per variant record as given in the original tsv file(s)
* `_all_vcf_only.tsv` - all variants present only in the vcfs (i.e. not present in the CombinedVariantOutput tsvs), this will have one line per sample per variant record as read from the vcf
* `_all_tsv_no_annotation.tsv` - all variants in the CombinedVariantOutput tsv files that have no proper annotation - this is crap left in by the local app that isn't used for interpretation, this is removed from the comparison of tsv vs vcf and output to a separate tsv