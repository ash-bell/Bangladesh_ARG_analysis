wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR02012.1.HMM

hmmsearch --tblout SAMPLE_recA.txt -E 1e-50 TIGR02012.1.HMM ../genecalls/SAMPLE_all_prot.faa;
# get genecall abundance from 04_genecalls_ARGs_mapping.sh line 44
