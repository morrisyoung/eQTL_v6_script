## build the qualified sample-tissue rep; to be further used later on
## two criterias: 1. sample has genotype information; and 2. sample is in eQTL (sample size>=60) tissue
## right now I use 60 as the threashold for eQTL tissues

## right after we get the sample list, we should process the rpkm file to only leave these qualified samples (for convenience of up-coming processing and computation)  --> done in this script


## notes:
##	1. first of all remove all samples of individuals that have no genotype information
##	2. then check how these samples distributed in all tissues, to decide how many tissues to involve
##	3. 





##=====================
##==== global variables
##=====================
individual_rep = {}		# hashing all the individuals with genotype information
sample_tissue_map = {}		# mapping all the samples into their tissue types
filter = 100			# TODO: we can change this to get more or less eTissues







##==================
##==== sub-routines
##==================
def get_individual_id(s):
	## naively find the second '-'
	id = ''
	count = 0
	for i in range(len(s)):
		if s[i] == '-':
			count += 1

		if count == 2:
			break

		id += s[i]

	return id





if __name__ == '__main__':



	# I need all the following (for the network model):
	'''
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_individuals_test
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_individuals_train
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count
	phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60
	'''






	##========================================================================================
	##==== get the sample-tissue mapping for all the samples
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type
	##========================================================================================
	"""
	file = open("../data_source/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	file1 = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'w')
	count = 0
	while 1:
		line = file.readline()[:-1]	# can't strip all "\t\t\t\t", as they are place holder
		count += 1
		if count <= 11:  ## 11
			continue

		if not line:
			break

		line = line.split('\t')
		sample = line[1]
		tissue1 = line[12]
		tissue2 = line[14]

		file1.write(sample + '\t' + tissue1 + '\t' + tissue2 + '\n')

	file.close()
	file1.close()
	"""






	##========================================================================================
	##==== remove samples that have no genotype information
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype
	##========================================================================================
	"""
	individuals = {}
	file = open("../data_source/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_sample_IDs.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		id = get_individual_id(line)
		individual_rep[id] = 1

	file.close()




	file = open("../data_source/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'w')
	file.readline()
	file.readline()

	# get the effective list for selection
	index_map = {}
	line = (file.readline()).strip()
	line = line.split('\t')
	file1.write(line[0] + '\t')
	for i in range(2, len(line)):
		sample = line[i]
		individual = get_individual_id(sample)
		if individual not in individual_rep:
			continue
		else:
			file1.write(line[i] + '\t')
			index_map[i] = 1
	file1.write('\n')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(2, len(line)):
			if i not in index_map:
				continue
			else:
				file1.write(line[i] + '\t')
		file1.write('\n')

	file.close()
	file1.close()
	"""



	



	##===================================================================================================
	##==== counting samples in each tissue
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_#size
	##===================================================================================================
	# sample_tissue_map
	"""
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'r')
	sample_tissue_map = {}
	while 1:
		line = file.readline()[:-1]	# can't strip '\t\t\t'
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_map[sample] = tissue
	file.close()

	# counting
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	print "there are",
	print len(sample_list),
	print "different samples from the rpkm file."

	counting = {}
	for sample in sample_list:
		tissue = sample_tissue_map[sample]
		if tissue not in counting:
			counting[tissue] = 1
		else:
			counting[tissue] += 1

	print "they are distributed among",
	print len(counting),
	print "different tissue types."

	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count", 'w')
	for tissue in counting:
		file.write(tissue + '\t' + str(counting[tissue]) + '\n')
	file.close()

	# filtering the counts
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_" + str(filter), 'w')
	count = 0
	for tissue in counting:
		if counting[tissue] >= filter:
			count += 1
			file.write(tissue + '\t' + str(counting[tissue]) + '\n')
	print "# of tissue with sample size >= " + str(filter) + ":",
	print count
	file.close()
	"""











	##======================================================================================================
	##==== extracting eQTL samples (the eTissue is defined as sample size >= #size)
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_#size_samples"
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample
	##======================================================================================================
	"""
	## get all the samples for tissue count >= filter
	## eQTL_tissue
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_" + str(filter), 'r')
	eQTL_tissue = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = (line.split('\t'))[0]
		eQTL_tissue[tissue] = []
	file.close()

	## sample_list
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	## sample_tissue_rep
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'r')
	sample_tissue_rep = {}
	while 1:
		line = file.readline()[:-1]
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_rep[sample] = tissue
	file.close()

	## fill in the eQTL_tissue rep
	for sample in sample_list:
		tissue = sample_tissue_rep[sample]
		if tissue in eQTL_tissue:
			eQTL_tissue[tissue].append(sample)

	# save the rep
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_" + str(filter) + "_samples", 'w')
	for tissue in eQTL_tissue:
		file.write(tissue + '\t')
		for sample in eQTL_tissue[tissue]:
			file.write(sample + '\t')
		file.write('\n')
	file.close()

	##============ process the rpkm matrix to get eQTL samples ==============
	## get the sample_rep first
	sample_rep = {}
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_" + str(filter) + "_samples", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')[1:]
		for sample in line:
			sample_rep[sample] = 1
	file.close()


	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'w')

	# filter all the samples again
	index_rep = {}
	line = (file.readline()).strip()
	line = line.split('\t')
	file1.write(line[0] + '\t')
	for i in range(1, len(line)):
		sample = line[i]
		if sample in sample_rep:
			index_rep[i] = 1
			file1.write(sample + '\t')
	file1.write('\n')
	
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(1, len(line)):
			if i in index_rep:
				file1.write(line[i] + '\t')
		file1.write('\n')

	file.close()
	file1.close()
	"""










	##======================================================================================================
	##==== remove all the NULL genes as defined (testing for all samples)
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize
	##======================================================================================================
	# GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'w')


	file.close()
	file1.close()






	##=============================================================================================
	##==== normalizing all the samples (here we use Log normalize other than the previous Quantile)
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize
	##=============================================================================================
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'w')


	file.close()
	file1.close()















	"""
	##======================================================================================================
	##==== separating all the esamples into their tissues
	##==== target: expression_by_etissue/tissue_list.txt
	##==== target: expression_by_etissue/tissue_x.txt
	##======================================================================================================
	# get the etissue list
	eQTL_tissue = {}
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_100_samples", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		eQTL_tissue[tissue] = {}
		for i in range(1, len(line)):
			sample = line[i]
			eQTL_tissue[tissue][sample] = 1

	# tissue list
	file = open("../data_processed/expression_by_etissue/tissue_list.txt", 'w')
	tissue_list = []
	count = 0
	for tissue in eQTL_tissue:
		tissue_list.append(tissue)
		count += 1
		file.write(tissue + '\t' + str(count) + '\n')
	file.close()

	# extract esamples into their etissues
	count = 0
	for tissue in tissue_list:
		count += 1
		file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'r')
		file1 = open("../data_processed/expression_by_etissue/tissue_" + str(count) + ".txt", 'w')

		# filter all the samples again
		index_rep = {}
		line = (file.readline()).strip()
		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(1, len(line)):
			sample = line[i]
			if sample in eQTL_tissue[tissue]:
				index_rep[i] = 1
				file1.write(sample + '\t')
		file1.write('\n')
		
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			file1.write(line[0] + '\t')
			for i in range(1, len(line)):
				if i in index_rep:
					file1.write(line[i] + '\t')
			file1.write('\n')

		file.close()
		file1.close()
	"""

