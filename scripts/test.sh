./scripts/ExcludonFinder \
  -f data/example/E.coli_K12_MG1655.fasta \
  -1 data/example/test_R1.fastq \
  -2 data/example/test_R2.fastq \
  -g data/example/E.coli_K12_MG1655.gff \
  -t 0.5 \
  -j 4 \
  -o test_output
