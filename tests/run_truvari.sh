#!/usr/bin/env bash

singularity exec docker://quay.io/biocontainers/truvari:5.3.0--pyhdfd78af_0 \
    truvari anno svinfo -o truvari.vcf test.vcf
