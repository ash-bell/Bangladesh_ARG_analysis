### contig taxnomy with kraken2 and kaiju ###

# get taxid of scaffolds
kaiju -z 16 \
 -t $WKDIR/ashley/tools/kaiju/kaijudb/nodes.dmp \
 -f $WKDIR/ashley/tools/kaiju/kaijudb/nr/kaiju_db_nr.fmi \
 -i $WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/assemblies/REPLACE_ME_scaffolds.fasta \
 -o REPLACE_ME_scaffolds_kaiju.out

# get taxid of scaffolds
kraken2 --db $WKDIR/ashley/tools/kraken2/2022_02_11_DB/kraken2_db/ \
--threads 16 \
--gzip-compressed \
$WKDIR/agb214/ARG_bangladesh_ponds/read_mapping/assemblies/REPLACE_ME_scaffolds.fasta > REPLACE_ME_scaffolds_kraken.out

# merge kaiju and kraken outputs using last common ancestor for ambigious reads
kaiju-mergeOutputs \
-i <(sort -k2,2 REPLACE_ME_scaffolds_kaiju.out) \
-j <(sort -k2,2 REPLACE_ME_scaffolds_kraken.out) \
-o REPLACE_ME_combined.out \
-c lca \
-t $WKDIR/ashley/tools/kaiju/kaijudb/nodes.dmp

# use taxonkit to convert NCBI taxIDs into taxonmic names
taxonkit lineage -c -i 3 -t REPLACE_ME_combined.out | awk -F'\t' 'BEGIN {OFS=FS} {gsub(-1,0,$4);print}' | taxonkit reformat -I 4 -t -F | cut -f 2,7 > REPLACE_ME_combined_rmft.tsv
