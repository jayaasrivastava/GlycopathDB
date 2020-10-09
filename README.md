Instructions to execute Annotation.py

Input files
	genome_filelist (a list of filenames for each protein sequences of genome obtained from NCBI)
	profile_info (profile_table.csv)
	seq_database (Database.fasta)
	pathway_master (Pathway_master.csv)
	reaction_master (Reaction_master.csv)
	pathway_sequence (Pathway_sequence.csv)
	blast_threshold (blast_thresholds.csv)
	concatenated_hmmfile (all_hmms)

Output files
	alert_output
	pathway_file
	pathway_proteins
	predictions_output

This code executes the following		       							   
	1)Obtains the file name of the genome from genome_filelist whose protein sequences are to be scanned 
	2)Opens the genome_file (cointaining protein sequences in fasta format) downloaded from NCBI
	3)Scans all profiles stored in concatenated_hmmfile against genome_file and stores results in hmmscan_output. Compares domain bit scores from hmmscan_output against threshold of each profile stored in profile_info. Hits are stored in hmm_dict
	4)Runs BLASTp to find homologs of sequences stored in seq_database from genome_file and compares sequence similarity and subject coverage that is obtained from blast_threshold. Stores hits in blast_dict
	5)Applies rules (described 10.1099/mgen.0.000452)
	6)Merges hmm_dict and blast_dict to merged_dict
	7)Checks if there are duplicates from HMM and BLAST (same protein a hit in both HMMscan and BLASTp scan) and writes to alert_output
	6)Scans pathway_master, pathway_sequence and reaction_master to find if homologs are present for each pathway in merged_dict. If yes, then it will publish that in pathway_file (pathways in a genome) and pathway_proteins (proteins associated with each step of a predicted pathway)             
	7)Additionally, it will write all predictions from HMM and BLASTp scan to predictions_output		       			   
