#!/bin/bash

mkdir data/genome/
cd data/genome/

for i in {1..20}; do
  curl -O "http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.0.dna.chromosome.${i}.fa.gz"
done

zcat Glycine_max.Glycine_max_v2.0.dna.chromosome.*.fa.gz > glymax.fa