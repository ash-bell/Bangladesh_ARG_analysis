# get taxid of forward and reverse reads
kaiju -z 16 \
 -t $WKDIR/ashley/tools/kaiju/kaijudb/nodes.dmp \
 -f $WKDIR/tools/kaiju/kaijudb/nr/kaiju_db_nr.fmi \
 -i $WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_qtrimmed.fwd.fq.gz \
 -j $WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_qtrimmed.rev.fq.gz \
 -o REPLACE_ME_qtrimmed_kaiju.out

# get taxid of unpaired reads
kaiju -z 16 \
 -t $WKDIR/ashley/tools/kaiju/kaijudb/nodes.dmp \
 -f $WKDIR/ashley/tools/kaiju/kaijudb/nr/kaiju_db_nr.fmi \
 -i $WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_merged.fq.gz \
 -o REPLACE_ME_merged_kaiju.out

# get taxid of forward and reverse reads
kraken2 --db $WKDIR/ashley/tools/kraken2/2022_02_11_DB/kraken2_db/ \
--threads 16 \
--gzip-compressed \
--paired \
$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_qtrimmed.fwd.fq.gz \
$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_qtrimmed.rev.fq.gz > REPLACE_ME_qtrimmed_kraken.out

# get taxid of unpaired reads
kraken2 --db $WKDIR/ashley/tools/kraken2/2022_02_11_DB/kraken2_db/ \
--threads 16 \
--gzip-compressed \
$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/read_taxonomy/REPLACE_ME_merged.fq.gz > REPLACE_ME_merged_kraken.out

# merge kaiju and kraken outputs using last common ancestor for ambigious reads
kaiju-mergeOutputs \
-i <(sort -k2,2 REPLACE_ME_qtrimmed_kaiju.out REPLACE_ME_merged_kaiju.out) \
-j <(sort -k2,2 REPLACE_ME_qtrimmed_kraken.out REPLACE_ME_merged_kraken.out) \
-o REPLACE_ME_combined.out \
-c lca \
-t $WKDIR/ashley/tools/kaiju/kaijudb/nodes.dmp

# use taxonkit to convert NCBI taxIDs into taxonmic names
taxonkit lineage -c -i 3 -t REPLACE_ME_combined.out | awk -F'\t' 'BEGIN {OFS=FS} {gsub(-1,0,$4);print}' | taxonkit reformat -I 4 -t -F | cut -f 2,7 > REPLACE_ME_combined_rmft.tsv
