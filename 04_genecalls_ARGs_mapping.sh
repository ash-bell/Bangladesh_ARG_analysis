### ARG search with RGI ###
rgi main \
-i assemblies/SAMPLE_scaffolds.fasta \
-o assemblies/SAMPLE_RGI_BLAST.txt \
-t contig \
-a BLAST \
-n 16 \
--low_quality \
-d NA \
--split_prodigal_job \
--exclude_nudge

### get gene abundance using bbmap ###
# map forward and reverse reads to gene calls from RGI
bbmap.sh ref=$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/genecalls/SAMPLE_nucl.fna \
build=GENECALLS \ # each gene called metagenome is referenced as a number (SAMPLE_nucl.fna). F1A = 1; F8D = 24. Allows indexing from the same folder
path=$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/protein_mappings/ \
in=$WKDIR/agb214/ARG_bangladesh_ponds/10218_SAMPLE/qtrimmed.fq.gz \
out=SAMPLE_genecalls_vs_SAMPLE_qtrimmed.bam

# map unpaired reads to gene calls from RGI
bbmap.sh ref=$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/genecalls/SAMPLE_nucl.fna \
build=GENECALLS \
path=$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/protein_mappings/ \
in=$WKDIR/agb214/ARG_bangladesh_ponds/10218_SAMPLE/merged.fq.gz \
out=SAMPLE_genecalls_vs_SAMPLE_merged.bam

#merge the unpaired and pair read abundances together into one bam file
samtools merge -o SAMPLE_genecalls_vs_SAMPLE_combined.bam \
SAMPLE_genecalls_vs_SAMPLE_merged.bam \
SAMPLE_genecalls_vs_SAMPLE_qtrimmed.bam

# sort the combined bam file
samtools sort -o SAMPLE_genecalls_vs_SAMPLE_combined.srt.bam \
SAMPLE_genecalls_vs_SAMPLE_combined.bam

# calculate gene coverage from mapping file using transcripts per kilobase million (TPM)
coverm contig \
--bam-files SAMPLE_genecalls_vs_SAMPLE_combined.srt.bam \
--min-covered-fraction 0.7 \ # this means 70% of a gene must be covered by reads to be included in downstream analysis
--min-read-percent-identity 0.95 \ # reads must have 95% with a gene call they are mapping against
--min-read-aligned-percent 0.9 \ # and 90% of the gene must map to a gene call
--methods tpm \
--threads 16 > SAMPLE_genecalls_vs_SAMPLE_combined.TPM.tsv
