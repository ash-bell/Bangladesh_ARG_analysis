#Interleave reads
reformat.sh in1=REPLACE_ME_r1.fq.gz in2=REPLACE_ME_r2.fq.gz out=interleave.fq.gz -Xmx1000g
ln -s interleave.fq.gz temp.fq.gz

#Remove optical duplicates
clumpify.sh in=temp.fq.gz out=clumped.fq.gz dedupe optical -Xmx1000g
rm temp.fq.gz; ln -s clumped.fq.gz temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=temp.fq.gz out=filtered_by_tile.fq.gz -Xmx1000g
rm temp.fq.gz; ln -s filtered_by_tile.fq.gz temp.fq.gz

#Trim adapters. Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=temp.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered threads=8 -Xmx1000g
rm temp.fq.gz; ln -s trimmed.fq.gz temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=temp.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality threads=8 -Xmx1000g
rm temp.fq.gz; ln -s filtered.fq.gz temp.fq.gz

#Decontamination by mapping can be done here.
#JGI removes these in two phases:
#1) common microbial contaminants (E.coli, Pseudomonas, Delftia, others)
#2) common animal contaminants (Human, cat, dog, mouse)

#Human
bbmap.sh build=1 path=$PATH/metag/ashley/eDNA/read_processing minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim in=temp.fq.gz outu=ex_human.fq.gz threads=8 -Xmx1000g
rm temp.fq.gz; ln -s ex_human.fq.gz temp.fq.gz

#Error-correct phase 1
bbmerge.sh in=temp.fq.gz out=ecco.fq.gz ecco mix vstrict ordered -Xmx1000g
rm temp.fq.gz; ln -s ecco.fq.gz temp.fq.gz

#Error-correct phase 2
clumpify.sh in=temp.fq.gz out=eccc.fq.gz ecc passes=4 reorder -eoom -Xmx1000g
rm temp.fq.gz; ln -s eccc.fq.gz temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
#Alternatively, bbcms.sh can be used if Tadpole still runs out of memory.
tadpole.sh in=temp.fq.gz out=ecct.fq.gz ecc k=62 ordered threads=8 -eoom -Xmx1000g
#bbcms.sh in=temp.fq.gz out=ecct.fq.gz ecc k=62 ordered threads=8 -eoom -Xmx1000g
rm temp.fq.gz; ln -s ecct.fq.gz temp.fq.gz

#Merge
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
bbmerge-auto.sh in=temp.fq.gz out=merged.fq.gz outu=unmerged.fq.gz strict k=93 extend2=80 rem ordered threads=8 -eoom -Xmx1000g

#Quality-trim the unmerged reads.
bbduk.sh in=unmerged.fq.gz out=qtrimmed.fq.gz qtrim=r trimq=10 minlen=70 ordered threads=8 -eoom -da -Xmx1000g

# --- Assembly ---
#Assemble with Spades
spades.py --meta -k 25,55,95,125 --only-assembler --phred-offset 33 -s merged.fq.gz --12 qtrimmed.fq.gz -o spades_out -m 1000 -t 8

#spades.py --restart-from last -o spades_out
