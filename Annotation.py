#/usr/local/bin/python

import subprocess
import sys
import os
import string
from sets import Set
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
import csv
import itertools

#dict element to store HMM and BLAST results. Initialized headers
checklist = []
hmm_dict = {}
blast_dict = {}
merged_dict={}	
profblast_res_dict = {} #this dict is only for all profiles and PUIDs which appear in results. It will be used for pathway prediction

def parse_and_check():

	genome_file_reader = csv.reader(genome_filelist, delimiter= ',')
	for file_row in genome_file_reader:
		genome_id = file_row[0]
		file_name = file_row[1]
		print(genome_id)
		genome_file = (genomefiles_path + file_name)
		run_hmmscan(genome_file,genome_id)
		run_blast(genome_file,genome_id)
		collate(genome_id)
		check_duplicates()
		predicted_pathways(genome_id)
	else:
		print("genome scanning successful")

####################################This function performs the following#####################################
#	1)Runs HMMSCAN command. It uses SearchIO parsing feature of biopython which is described in 	    	#
#	  section 1. All profiles which appear in HMMscan are verified for their threshold in   	    		#	
#	  compare_threshold(). If satisfied, it will fill two dictionaries:				    					#
#		1)prof_blast_dict: which only contains genome_id:acc_id as key and hit profiles as values   		#
#		2)hmm_dict: which contains, in addition to above key and hit,aln range and bit score of hit 		#
#############################################################################################################

def run_hmmscan(genome_file,genome_id):

        command_hmm = ["hmmscan", "--domtblout", hmmscan_output, "--noali", concatenated_hmmfile, genome_file]
        FNULL = open(os.devnull, 'w')
        subprocess.call(command_hmm, stdout=FNULL, stderr=subprocess.STDOUT)

        with open(hmmscan_output, 'rU') as input:
                num_hits = 0
                for qresult in SearchIO.parse(input, 'hmmscan3-domtab'):
                        hits = qresult.hits
                        num_hits = len(hits)
                        if num_hits > 0:
                                for hsp in hits:
                                        for frag in hsp:
						profile_id = hsp.id
                                                identification = hsp.query_id
                                                accession = identification.split("|")[1] if "|" in identification else (identification)  #accomodate differing seq headers in uniprot and other websites
						bit_score = int(frag.bitscore)
                                                q_start = str(frag.query_start + 1)                      		        
                                                q_end = str(frag.query_end)								 #in hmmsearch, this is hit_start/end			        
                                                q_range = q_start + "-" + q_end
						profile_item_values = profile_id + "," + q_range + "," + str(bit_score)
						hmm_dict_key = genome_id + "_" + accession
						if compare_threshold(profile_id,bit_score) == 'Y':
							profblast_res_dict.setdefault(hmm_dict_key, [])			#D110419 (This step will append all matching profiles to a protein)
	                                                hmm_dict.setdefault(hmm_dict_key, [])
	                                                hmm_info_list = []
	                                                hmm_info_list = [profile_id,bit_score,q_range]
	                                                hmm_dict[hmm_dict_key].append(hmm_info_list)
							profblast_res_dict[hmm_dict_key].append(profile_id)
						else:
							pass
	print("hmm complete")


def compare_threshold(profile_id,bit_score):

	profile_table_reader = csv.reader(profile_info, delimiter= ',')
       	for row in profile_table_reader:
		if profile_id in row:
			threshold = int(row[4])
			break
	if (bit_score >= threshold):
		return 'Y'
	else:
		return 'N'
				

####################################This function performs the following#####################################
#	1)Searches for homologous entries in sequence database								 				    #	
#	2)A list of PUIDs and their cut off ssim and scov values are stored in blast_thresholds.csv			 	#
#	2)BLAST results will be parsed and if the target id to which the genome accession found homology	    #
#	  is present in above file, the correponding thresholds will be matched and written in blast_dict   	#
#	5)Two master dictionaries are used below															    #
#	   1)blast_dict which will contain output with aln_range, ssim and scov values							#
#	   2)prof_blast_dict(also populated in HMM scan) with only hit PUIDs									#
#############################################################################################################

def run_blast(genome_file,genome_id):

        command_blast = ["blastp", "-db", seq_database, "-word_size", "2", "-query", genome_file, "-out", blast_output, "-outfmt", "6 qseqid sseqid ppos slen qstart qend sstart send"]
        subprocess.call(command_blast)

	threshold_puids = []
	btreader = csv.reader(blast_threshold, delimiter= ',')
	for row1 in btreader:
		threshold_puids.append(row1[0])	
	
        if os.path.getsize(blast_output) > 0:
                with open(blast_output) as infile:
                        for line in infile:
                                entry = line.split("\n")[0]                                                   #as output from blastp has an extra newline character in the end
                                qid = entry.split("\t")[0]
                                tid = entry.split("\t")[1]
                                ssim = float(line.split("\t")[2])
                                subject_length = float(entry.split("\t")[3])
                                subject_start = float(entry.split("\t")[6])
				subject_end = float(entry.split("\t")[7])
				qstart = entry.split("\t")[4]
				qend = entry.split("\t")[5]
				aln_range = qstart + "-" + qend
                                subject_cov = float((subject_end - subject_start + 1)/subject_length) * 100
				if subject_cov >= 100:
					 subject_cov = 100
				else:
					subject_cov = int(subject_cov)
				blast_dict_key = genome_id + "_" + qid
				if tid in threshold_puids:
					btreader = csv.reader(blast_threshold, delimiter= ',')
					for row1 in btreader:
						if tid in row1:
							blast_ssim_threshold = (row1[1])
							blast_scov_threshold = (row1[2])
							if (ssim >= float(blast_ssim_threshold)) and (int(subject_cov) >= int(blast_scov_threshold)): 
								profblast_res_dict.setdefault(blast_dict_key, [])					
								blast_dict.setdefault(blast_dict_key, [])
								profblast_res_dict[blast_dict_key].append(tid)							
					       			blast_desc_dict = {}
								blast_info_list = []
								blast_info_list = [tid,ssim,aln_range,subject_cov]
					     			blast_dict[blast_dict_key].append(blast_info_list)
								blast_desc_dict = {}
				else:
					if subject_cov >= 90 and ssim >= 50:
						blast_dict.setdefault(blast_dict_key, [])
						profblast_res_dict.setdefault(blast_dict_key, [])
						profblast_res_dict[blast_dict_key].append(tid)					
			       			blast_desc_dict = {}
						blast_info_list = []
						blast_info_list = [tid,ssim,aln_range,subject_cov]
			     			blast_dict[blast_dict_key].append(blast_info_list)
						blast_desc_dict = {}

	print("blast ended")	

####################################This function performs the following#########################################
#					RULES!!! 																					#
#	1)If any of the rules (described below) is not met, the exception will be written in alert_output			#
#	2)Post the rules, prof_blast_dict and hmm_dict will be read. If the length of values (number of 			#
#	  hit_profiles is more than one in prof_blast_dict), (for eg. aminotransferases), it will pick the most 	#
#	  specific one by verifying it from substrate_specific_hmms list. If its present in this list, then it  	#
#	  take the corresponding values from hmm_dict for that accession. 											#
#	3)It will perform a check if its a case of fusion protein by comparing the start and end amino acid     	#
#	   (Note: this is done by first making a start and end list, and then checking if ending aa is less			#
#	   than starting amino acid of next domain)																	#
#	   A list is made for comparison anticipating cases of more than two fusion domains. For instance			#
#	   xyz[305-500], xyz[3-300], xyz[501-800] are picked by three different profiles.							#
#	   A list will be created with start and end amino acids, such that the two lists are 						#
#			start_range_list = [305,3,501]																		#
#			end_range_list = [500,300,800]																		#
#	   The list will be sorted and then running a loop 3-1 times(len of list - 1), it will check if 			#
#	   end_range_list[0] is less than start_range_list[1], if yes, then it will print in merged_dict			#
#       4)Note: junk values 1111, 2222 and 3333 are filled in NA integer columns that are replaced by blanks in #
#         output files generated                                                                                #
#################################################################################################################

def collate(genome_id):

	specific_aminoT_hmms = ['GPE01430','GPE01710','GPE01530']
	C4_C3_aminoTs = ['GPE01910','GPE01830']
	uridylylT = ['GPE00430','GPE00530']
	guanylylT = ['GPE00620','GPE00720']
	mutases = ['GPE09630','GPE09330']
	isomerases = ['GPE07030', 'GPE07130', 'GPE07230', 'GPE07330']
	dehydratase1 = ['GPE05332', 'Q81A42_1-328']				#used in check_duplicates()
	dehydratase2 = ['GPE05331', 'GPE05332']
	dehydratase3 = ['GPE05430', 'GPE05510']
	dehydratase4 = ['GPE05430', 'GPE02230']
	dehydrogenases1 = ['GPE03210', 'GPE03130']
	dehydrogenases2 = ['GPE03030', 'GPE03430']
	spc_mutases = ['GPE09130', 'GPE09230']
	c2_epimerases = ['GPE02030', 'GPE02110']
	

	blasthmm_hits = list(profblast_res_dict.values())
	flat_list = [y for x in blasthmm_hits for y in x]

	for hmm_accession,profile_item in hmm_dict.items():
		value_len = len([item for item in profile_item if item])
		if value_len > 1:
			domain_range_list = []
			for hmm_item in profile_item:
				start = int(hmm_item[2].split("-")[0])
				end = int(hmm_item[2].split("-")[1])
				aligned_range = set(range(start,end))	
				domain_range_list.append(aligned_range)
			val = profblast_res_dict.get(hmm_accession)
			collated_domains=list(map(list, {tuple(i for i, ss in enumerate(domain_range_list) if set(s).intersection(ss)) for s in domain_range_list}))
			for items_indices in collated_domains:
				prof_domains = []
				for item_index in items_indices:
					dom_profile_item = val[item_index]
					prof_domains.append(dom_profile_item)
				#rules
				common1 = [i for i in prof_domains if i in specific_aminoT_hmms]	
				common2 = [i for i in prof_domains if i in C4_C3_aminoTs]
				common3 = [i for i in prof_domains if i in uridylylT]
				common4 = [i for i in prof_domains if i in guanylylT]
				common5 = [i for i in prof_domains if i in mutases]
				common6 = [i for i in prof_domains if i in isomerases]
				common7 = [i for i in prof_domains if i in dehydratase2]
				common8 = [i for i in prof_domains if i in dehydratase3]
				common9 = [i for i in prof_domains if i in dehydratase4]
				common10 = [i for i in prof_domains if i in dehydrogenases1]
				common11 = [i for i in prof_domains if i in dehydrogenases2]
				common12 = [i for i in prof_domains if i in spc_mutases]
				common13 = [i for i in prof_domains if i in c2_epimerases]	

				if len(common1) > 0:
					if len(common2) > 0:
						if 'GPE01230' in prof_domains:
							pass
						else:
							message = (str(common2) + ' found but not GPE01230')
							write_output(hmm_accession,message)
					else:
						message = (str(common1) + ' found but not GPE01830/GPE01910')
						write_output(hmm_accession,message)
				if len(common2) > 0:
					if 'GPE01230' in prof_domains:
						pass
					else:
						message = (str(common2) + ' found but not GPE01230')
						write_output(hmm_accession,message)
				if 'GPE01430' in prof_domains:
					idx = val.index('GPE01430')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE01530' in prof_domains:
					idx = val.index('GPE01530')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE01710' in prof_domains:
					idx = val.index('GPE01710')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE01830' in prof_domains and len(common1) == 0:	
					idx = val.index('GPE01830')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE01910' in prof_domains and len(common1) == 0:
					idx = val.index('GPE01910')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE00210' in prof_domains:
					if 'GPE00231' in prof_domains:
						pass
					else:
						message = ('GPE00210 found but not GPE00231')
						write_output(hmm_accession,message)
					idx = val.index('GPE00210')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE02430' in prof_domains:
					if 'GPE02530' in prof_domains:
						pass
					else:
						message = ('GPE02430 found but not GPE02530')
						write_output(hmm_accession,message)
					idx = val.index('GPE02430')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE03210' in prof_domains:
					if 'GPE03430' in prof_domains:
						pass
					else:
						message = ('GPE03210 found but not GPE03430')
						write_output(hmm_accession,message)
					idx = val.index('GPE03210')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE03130' in prof_domains:
					if 'GPE03430' in prof_domains:
						pass
					else:
						message = ('GPE03130 found but not GPE03430')
						write_output(hmm_accession,message)
					idx = val.index('GPE03130')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE09130' in prof_domains:
					if 'GPE09630' in prof_domains:
						pass
					else:
						message = ('GPE9130 found but not GPE09630')
						write_output(hmm_accession,message)
					idx = val.index('GPE09130')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif 'GPE09230' in prof_domains:
					if 'GPE09630' in prof_domains:
						pass
					else:
						message = ('GPE9230 found but not GPE09630')
						write_output(hmm_accession,message)
					idx = val.index('GPE09230')
					hmm_list = profile_item[idx]
					merged_dict_create(hmm_accession,hmm_list[0],hmm_list[1],hmm_list[2],"NA",1111,2222)
				elif len(common2) > 1:
					if len(set(prof_domains)) == len(prof_domains):		#added on D19/11/19 here as well as in subsequent rules to deal with cases where single ORF can have identical domains fused. For instance, if same protein has two domains hit to GPE01910, then previous condition (len(common2) > 1)) is incomplete. Nonetheless, since both the domains are to be written in output, thw write para i.e. (for hmm_item in profile_item:) was outside this condition 
						message = ('Hit for both GPE01910 and GPE01830')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common3) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE00430 and GPE00530')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common4) > 1:
					if len(set(prof_domains)) == len(prof_domains):	
						message = ('Hit for both GPE00620 and GPE00720')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common5) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE09630 and GPE09330')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common6) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for more than one isomerase profiles')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common7) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE05331 and GPE05332')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common8) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE05510 and GPE05430')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common9) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE05430 and GPE02230')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common10) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE03210 and GPE03130')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common11) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE03030 and GPE03430')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common12) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE09130 and GPE09230')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(common13) > 1:
					if len(set(prof_domains)) == len(prof_domains):
						message = ('Hit for both GPE02030 and GPE02110')
						write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				elif len(prof_domains) > 1:										# for other cases which have multiple hits for one domain not
					message = ('Hit for more than one profiles or has multiple domains')				#covered in aforementioned rules
					write_output(hmm_accession,message)
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
				else:
					for hmm_item in profile_item:
						merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)
		else:
			for hmm_item in profile_item:
				merged_dict_create(hmm_accession,hmm_item[0],hmm_item[1],hmm_item[2],"NA",1111,2222)

	enzyme_set1 = ['Q9XC60_1-372','P95700_1-371']
	enzyme_set2 = ['Q66DP5_1-329','P26395_1-330']
	enzyme_set3 = ['Q0P8T8_1-323','G8UHQ7_1-325']
	enzyme_set4 = ['O52793_1-471','Q9ALN6_1-486']
	enzyme_set5 = ['A9IH93_1-190','G3XD01_1-191']
	enzyme_set6 = ['O35826_410-684','Q9Y223_406-722']
	enzyme_set7 = ['O25093_209-517','A0A1D3UJQ1_1-321','Q0P8U5_1-274','Q3ESA1_1-366']
	enzyme_set8 = ['P71063_1-216','B0V6J7_1-219','Q0P9D1_1-195','Q5FAE1_197-403']

	for blast_accession,blast_puid in blast_dict.items():
		value_len = len([item for item in blast_puid if item])
		if value_len > 1:
			domain_range_list = []
			for blast_item in blast_puid:
				start = int(blast_item[2].split("-")[0])
				end = int(blast_item[2].split("-")[1])
				aligned_range = set(range(start,end))	
				domain_range_list.append(aligned_range)
			collated_domains=list(map(list, {tuple(i for i, ss in enumerate(domain_range_list) if set(s).intersection(ss)) for s in domain_range_list}))
			val = profblast_res_dict.get(blast_accession)
			for items_indices in collated_domains:
				blast_domains = []
				for item_index in items_indices:
					dom_blast_item = val[item_index]
					blast_domains.append(dom_blast_item)
				common1 = [i for i in blast_domains if i in enzyme_set7]
				common2 = [i for i in blast_domains if i in enzyme_set8]
				if(all(x in enzyme_set1 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif(all(x in enzyme_set2 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif(all(x in enzyme_set3 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif(all(x in enzyme_set4 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif(all(x in enzyme_set5 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif(all(x in enzyme_set6 for x in blast_domains)):
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif len(common1) > 1:
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				elif len(common2) > 1:
						blast_item = blast_puid[0]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
				else:
					for item_index in items_indices:
						blast_item = blast_puid[item_index]
						merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])
		else:	
			for blast_item in blast_puid:
				merged_dict_create(blast_accession,"NA",3333,blast_item[2],blast_item[0],blast_item[3],blast_item[1])

	pdeg = 'Q81A42_1-328'
	preq = 'Q81A43_1-303'
	if pdeg in flat_list and preq not in flat_list:
		message = 'pdeg found but not preq'
		write_output(genome_id,message)
	elif preq in flat_list and pdeg not in flat_list:
		message = 'preq found but not pdeg'
		write_output(genome_id,message)
	
	print("merged_dict created")
			
#########################################################################################################################
#	Merged_dict is the final dictionary that will also be output in Prediction_output file				#
#	The key is created by adding range (as we do in PUIDs) as key and then profile_hit, bit score, PUID, ssim,qcov  #
#	as values. 													#
#########################################################################################################################
						
def merged_dict_create(key,hit_profile,bit_score,aln_range,loner_id,scov,sim):

	new_key = key + ":" + aln_range
	item_str = [hit_profile,str(bit_score),loner_id,str(scov),str(sim)]
	if new_key in merged_dict.keys():
		pass
	else:
		merged_dict.setdefault(new_key, [])
		for i in item_str:
			merged_dict[new_key].append(i)

def write_output(reference,message):

       	with open(alert_output, 'a') as outfile:
	        outfile.write(reference + ":" + "\t" + message + '\n')
		
def check_duplicates():

	dehydratase1 = ['GPE05332', 'Q81A42_1-328']
	for key1,value1 in hmm_dict.items():
		if key1 in blast_dict:
			start_range_list = []
			end_range_list = []
			for key2,value2 in merged_dict.items():
				key_acc = key2.split(":")[0]
				key_range = key2.split(":")[1]
				if key1 == key_acc:
					start_range_list.append(int(key_range.split("-")[0]))
					end_range_list.append(int(key_range.split("-")[1]))
				hits = set(value2)
			start_range_list.sort()
			end_range_list.sort()
			for i in range(len(start_range_list)-1):
				if (end_range_list[i]) <= (start_range_list[i+1]):
					pass
				else:
					if set(dehydratase1) == hits:
						message = ('Hit to both GPE05332 and Q81A42_1-328')
						write_output(key1,message)
					message = 'found in both blast and HMM'
					write_output(key1,message)
					break


##added include pathway prediction feature (searchfor D110419)####################

def predicted_pathways(genome_id):	

	print("pathway processing started")			

	pw_master_reader = csv.reader(pathway_master, delimiter= ',')
	for pmrow in pw_master_reader:
		pwid = pmrow[0]								#storing pwid which will be used to fetch stepid, pw_seq_id and steptarget
		stepid_list = []
		pw_seq_list = []							#storing all step ids of a pathway, reinitialized here everytime a new pwid is read
		steptarget_list = []							#storing all step names of a pathway, reinitialized here everytime a new pwid is read	
		pw_sequence_reader = csv.reader(pathway_sequence, delimiter= ',')				
		for psrow in pw_sequence_reader:
			if psrow[1] == pwid:
				stepid_list.append(psrow[3])				#this will be used to find matching profile/BLAST id from reaction_master in line 320, 324 and 326
				pw_seq_list.append(psrow[0])				#this will be used to display pw_seq in pathway_proteins output
		for s in stepid_list:
			reaction_master_reader = csv.reader(reaction_master, delimiter= ',')
			for smrow in reaction_master_reader:
				if s == smrow[0]:
					if not smrow[2]:						#check if Profile hit colm is empty. Populate blast-id if yes
						if '/' in smrow[3]:					#checks if there are multiple BLAST hit ids in step master
							smrow_blast = smrow[3].split("/")								
							steptarget_list.append(smrow_blast)		#create a nested list of multiple profile/BLAST hit ids inside list steptarget_list
						else:
							steptarget_list.append(smrow[3])
					elif not smrow[3]:						#check if blast column is empty. Populate profile-id if yes
						if '/' in smrow[2]:					#checks if there are multiple profiles hit ids in step master
							smrow_profiles = smrow[2].split("/")
							steptarget_list.append(smrow_profiles)		#create a nested list of multiple profile/BLAST hit ids inside list steptarget_list	
						else:
							steptarget_list.append(smrow[2])
					else:								#This is true if there are hits for both prof and blast. for eg. Fructose-6-phosphate isomerase
						temp_hit_list = []					
						if '/' in smrow[3]:				
							smrow_blast = smrow[3].split("/")								
							temp_hit_list.append(smrow_blast)
							flat_list = [item for sub_list in temp_hit_list for item in sub_list]
						else:
							temp_hit_list.append(smrow[3])
							flat_list = temp_hit_list
						if '/' in smrow[2]:
							smrow_profiles = smrow[2].split("/")
							temp_hit_list.append(smrow_profiles)
							flat_list = [item for sub_list in temp_hit_list for item in sub_list]
						else:
							temp_hit_list.append(smrow[2])
							flat_list = temp_hit_list
						steptarget_list.append(flat_list)		

###################################################################################################
#following will create a reverse dictionary of prof_blast_dict. Basically one to many relationship#
#between accession and profile/PUID hit is reversed to one to many relationship between profile/  #
#PUID and accession. This will allow all hits for a particular pathway step to be obtained.       #
###################################################################################################
		accession_hit_dict = {}		
		for query_acc,hit in merged_dict.items():
			if hit[0] == 'NA':
				accession_hit_dict.setdefault(hit[2], [])
				accession_hit_dict[hit[2]].append(query_acc[16:])			
			elif hit[2] == 'NA':
				accession_hit_dict.setdefault(hit[0], [])
				accession_hit_dict[hit[0]].append(query_acc[16:])			
	
		pathway_accession = {}
		for x in steptarget_list:
			is_present = False
			if isinstance(x, list):
				for y in x:
					if y in accession_hit_dict:
						idx = steptarget_list.index(x)
						if idx in pathway_accession.keys():
							pass
						else: 
							pathway_accession[idx] = {}
						accession_id = accession_hit_dict.get(y)
						pathway_accession[idx][y] = accession_id
						is_present = True
			else:
				if x in accession_hit_dict: 
					idx = steptarget_list.index(x)
					pathway_accession[idx] = {}
					accession_id = accession_hit_dict.get(x)
					pathway_accession[idx][x] = accession_id
					is_present = True

		pathway_file = (genome_id + "_predicted_pathways")
		pathway_proteins = (genome_id + "_pathway_proteins")

		if len(pathway_accession) == len(steptarget_list):						
			with open(pathway_file, 'a+') as genome_pathway_file:
				genome_pathway_file.write(genome_id + ',' + pwid + '\n')
			with open(pathway_proteins, 'a+') as pathway_proteins_file:
				for j in range(len(steptarget_list)):					#this logic starts with looping over steptarget_list since it contains sequential profile/BLAST
					for step, res in pathway_accession.items():				# id as found in the pathway	
						hit_list = []
						hit_acc_list = []								
						for hit, hit_acc in res.items():
							hit_list.append(hit)
							hit_acc_list.append(hit_acc)	
						flat_list = [item for sublist in hit_acc_list for item in sublist]
						if(all(x in steptarget_list[j] for x in hit_list)):				
							pathway_proteins_file.write(genome_id + ',' + str(flat_list) + ',' + pw_seq_list[j] + '\n')

	prediction_output = (genome_id + "_genome_predictions")

	with open(prediction_output, 'a+') as outfile:
		for key,value in merged_dict.items():
			outfile.write(key + ',' + str(value) + '\n')

	profblast_res_dict.clear()
	hmm_dict.clear()
	blast_dict.clear()
	merged_dict.clear()

	print("pathway processing ended")


#########################-----------MAIN-------------############################
#################################################################################   

parse_and_check()