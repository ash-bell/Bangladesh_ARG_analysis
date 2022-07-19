### reassemble metagenome with metaplasmidspades for plasmids ###
metaplasmidspades.py \
--only-assembler \
--phred-offset 33 \
-s $WKDIR/ashley/10218_SAMPLE/merged.fq.gz \
--12 $WKDIR/ashley/10218_SAMPLE/qtrimmed.fq.gz 
-o $WKDIR/ashley/10218_SAMPLE/metaPlasmidSPAdes_output \
-m 1000 \
-t 8

#metaplasmidspades.py --continue -o $WKDIR/ashley/10218_SAMPLE/metaPlasmidSPAdes_output

### detect plasmids from assembly using PlasClass ###
python classify_fasta.py \
-f $WKDIR/ashley/plasmids/SAMPLE_plasmids_scaffolds.fasta \
-o SAMPLE_plasmids.list \ # keep any contig with score > 0.5
-p 16

### use virsorter pipeline to detect viruses and remove from downstream analysis ### 
virsorter run --keep-original-seq \
-i SAMPLE_plasmids_plasclass.fasta \
-w SAMPLE_VS2/vs2-pass1 \
--include-groups dsDNAphage,ssDNA \
--min-length 5000 \
--min-score 0.5 \
-j 16 \
all

checkv end_to_end \
SAMPLE_VS2/vs2-pass1/final-viral-combined.fa \
SAMPLE_VS2/checkv \
-t 16 \
-d $WKDIR/ashley/tools/checkv-db-v1.0

# remove any contig in "keep 1" / "keep 2" or "manual check" catagory. 
# See here https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=4

### detect ARGs in plasmids using RGI ###
rgi main \
-i SAMPLE_plasmids_not_viral.fasta \
-o SAMPLE_plasmids_ARGs.txt \
-t contig \
-a BLAST \
-n 16 \
--clean \
--low_quality \
-d plasmid \
--split_prodigal_job \
--exclude_nudge

### get identity of genecalls on plasmids
diamond blastx \
--query interesting_plasmid.fnn \
--db $WKDIR/dbs/nr/nr.2021-07-09.diamond.2.0.11.dmnd \
--out dmd_plasmid.out \
--outfmt 6 \
--threads 16
