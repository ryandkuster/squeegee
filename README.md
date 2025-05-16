# squeegee (for vg)

Squeegee is currently a rudimentary script that takes as input a vg deconstruct vcf file that has been normalized using bcftools norm.

The output is largely meant to allow users to have an INFO field annotated with SVTYPE, SVLEN, and END, as is common in structural variant calling tools (e.g., svim, sniffles, cutesv).

The currently supported SVTYPES are SNP, MNP, INS, and DEL. If you have an example vcf that has DUP/INV/BND, please send it my way and I'll attempt to add this ability if possible (rkuster@utk.edu).

## Input file

The input vcf must first be normalized using bcftools. See this [github documentation](https://github.com/noecochetel/North_American_Vitis_Pangenome/blob/main/0.07_infer_variants.md) and the corresponding [Cochetal et al. 2023 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03133-2).

```bash
bcftools norm --output vg_deconstruct_norm.vcf.gz --output-type v -m-any --fasta-ref ref.fna --check-ref w vg_deconstruct.vcf.gz
```

## Running squeegee

To run from the squeegee with default settings (no filtering, just add INFO fields):

```bash
./squeegee -v vg_deconstruct_norm.vcf.gz -o vg_squeegee.vcf.gz
```

## Filtering with squeegee

Similar to vcfbub, filter variants based on snarl level (LV) or other attributes.

`-m` is the minimum difference (SVLEN) between ALT and REF sequence lengths to keep.
`-M` is the maximum difference (SVLEN) between ALT and REF sequence lengths to keep.
`-t` is are the types of variants to keep (default is all).
`-l` is the maximum snarl level to keep.

```bash
./squeegee --help
```

...yields:

```
optional arguments:
  -h, --help    show this help message and exit
  -v V          Path to bcftools norm biallelic vg deconstruct file.
  -o O          Path to output vcf.
  -t T [T ...]  Space separated list of variant types to keep.
  -m M          Minimum SVLEN to keep.
  -M M          Minimum SVLEN to keep.
  -l L          Maximum level in snarl tree (LV) to keep.
```
