## this is used to initialzie the cis- parameters (actually learn these parameters from multi-linear regression)
## the learned parameters have no tissue specificity (that's why this is only the initialization)
## the results are saved in "../para_init/para_init_train_cis.txt", each line is for one gene (several SNPs are in the cis- region of this gene)
## we can later on use the same script to test the prediction precision, and draw the plot for all the ~20,000 genes


import time
import numpy as np


# global variables definition and initialization
num_gene = 0			# TBD
num_individual = 0		# TBD

# genotype:
snp_pos_list = []		# the position of all SNPs
snp_dosage_rep = {}		# load all the dosage (for cluster)

# expression:
sample_rep = {}			# hashing all training samples into their rpkm array
sample_list = []		# in-order sample list
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)

# information table:
gene_tss = {}			# TSS for all genes (including those pruned genes)
gene_xymt_rep = {}		# map all the X, Y, MT genes
gene_cis_index = {}		# mapping the gene to cis snp indices (start position and end position in the snp vector)

# result table:
para_rep = {}




def snp_dosage_load():
	global snp_dosage_rep

	for individual in snp_dosage_rep:
		for i in range(22):
			chr = str(i+1)
			snp_dosage_rep[individual].append([])
			file = open("../../GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/genotype_imputed/genotype_450_dosage_matrix_qc/chr" + chr + "/SNP_dosage_" + individual + ".txt", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break
 
				dosage = float(line)
				snp_dosage_rep[individual][i].append(dosage)
			file.close()
	return




if __name__ == "__main__":

	print "working on the multi-linear regression of cis- SNPs for all genes"
	time_start = time.time()


	##===================================================== genotype =====================================================
	print "snp_dosage_rep..."
	file = open("../data_processed/list_individuals.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line
		snp_dosage_rep[individual] = []
	file.close()

	num_individual = len(snp_dosage_rep)
	snp_dosage_load()


	## snp_pos_list
	print "snp_pos_list..."
	snp_pos_list = []
	for i in range(22):
		chr = str(i+1)
		snp_pos_list.append([])
		file = open("../../GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/genotype_imputed/genotype_450_dosage_matrix_qc/chr" + chr + "/SNP_info.txt", 'r')

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			snp = line[0]
			pos = int(line[1])
			snp_pos_list[i].append(pos)
		file.close()


	##===================================================== expression =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get training samples..."
	file = open("../data_processed/list_samples_train.txt", 'r')
	sample_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]


		for i in range(len(sample_list)):
			sample = sample_list[i]
			sample_rep[sample] = []

	file.close()


	##=================== query the RPKM according to the above list
	print "get expression matrix for these training samples..."
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')

	###
	file.readline()
	file.readline()
	sample_list = ((file.readline()).strip()).split('\t')[2:]
	index_rep = {}
	for i in range(len(sample_list)):
		sample = sample_list[i]
		if sample in sample_rep:
			index_rep[i] = sample

	###
	gene_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		gene_list.append(gene)
		rpkm_list = map(lambda x: float(x), line[2:])

		for i in range(len(rpkm_list)):
			rpkm = rpkm_list[i]
			if i in index_rep:
				sample = index_rep[i]
				sample_rep[sample].append(rpkm)
	file.close()

	###
	gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
	for i in range(len(gene_list)):
		gene = gene_list[i]
		gene_index_map[gene] = i

	###
	num_gene = len(gene_list)

	###
	sample_list = []
	for sample in sample_rep:
		sample_list.append(sample)

	###
	expression_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		expression_matrix.append(sample_rep[sample])
	expression_matrix = np.array(expression_matrix)

	###
	gene_tss = {}
	file = open("../data_processed/gene_tss.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		tss = int(line[2])
		gene_tss[gene] = (chr, tss)
	file.close()

	###
	gene_xymt_rep = {}
	file = open("../data_processed/gene_xymt.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		gene = line
		gene_xymt_rep[gene] = 1
	file.close()



	##===================================================== cis- region definition =====================================================
	# gene_cis_index
	for i in range(len(gene_list)):
		gene = gene_list[i]
		if gene in gene_xymt_rep:
			continue
		chr = int(gene_tss[gene][0])
		tss = gene_tss[gene][1]
		flag1 = 0
		flag2 = 0
		start = 0
		end = 0
		for j in range(len(snp_pos_list[chr-1])):
			if flag1 == 0:
				if (snp_pos_list[chr-1][j] - tss >= -1000000) and (snp_pos_list[chr-1][j] - tss <= 1000000):
					start = j
					flag1 = 1
			if flag1 == 1 and flag2 == 0:
				if (snp_pos_list[chr-1][j] - tss >= -1000000) and (snp_pos_list[chr-1][j] - tss <= 1000000):
					end = j
				else:
					flag2 = 1;
			if (flag1 == 1) and (flag2 == 1):
				break
		gene_cis_index[gene] = (start, end)
	file = open("../data_processed/gene_cis_range.txt", 'w')
	for gene in gene_cis_index:
		start = gene_cis_index[gene][0]
		end = gene_cis_index[gene][1]
		file.write(gene + '\t' + str(start) + '\t' + str(end) + '\n')
	file.close()



	gene_cis_index = {}
	file = open("../data_processed/gene_cis_range.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		start = int(line[1])
		end = int(line[2])
		gene_cis_index[gene] = (start, end)
	file.close()



	##===================================================== regression across all genes =====================================================
	for i in range(len(gene_list)):

		gene = gene_list[i]
		print gene

		if gene in gene_xymt_rep:
			continue

		chr = int(gene_tss[gene][0])
		start = gene_cis_index[gene][0]
		end = gene_cis_index[gene][1]

		##
		expression_array = []
		for j in range(len(sample_list)):
			rpkm = expression_matrix[j][i]
			expression_array.append(rpkm)
		expression_array = np.array(expression_array)

		##
		genotype_matrix = []
		for j in range(len(sample_list)):
			sample = sample_list[j]
			genotype_matrix.append([])
			individual = sample[:9]
			for k in range(start, end+1):
				dosage = snp_dosage_rep[individual][chr-1][k]
				genotype_matrix[j].append(dosage)
			genotype_matrix[j].append(1)  # we need the intercept
		genotype_matrix = np.array(genotype_matrix)
		## sample:
		#X = np.array([[1,2,3,1], [2,4,6,1], [3,6,9,1]])  # 1x1 + 2x2 + 3x3 + 1
		#y = np.array([14, 28, 42])
		#y = np.array([15, 29, 44])
		#m = np.linalg.lstsq(X, y)[0]
		try:
			m = np.linalg.lstsq(genotype_matrix, expression_array)[0]
			para_rep[gene] = m  ## there is an extra intercept here!!!
		except ValueError:
			print genotype_matrix
			## write the matrix into a file
			file = open("./temp/" + gene + ".txt", 'w')
			for i in range(len(genotype_matrix)):
				for j in range(len(genotype_matrix[i])):
					dosage = genotype_matrix[i][j]
					file.write(str(dosage) + '\t')
				file.write('\n')
			file.close()


	##===================================================== save all learned parameters =====================================================
	file = open("../para_init/para_init_train_cis.txt", 'w')
	for gene in para_rep:
		para_list = para_rep[gene]
		file.write(gene + '\t')
		for i in range(len(para_list)):
			para = para_list[i]
			file.write(str(para) + '\t')
		file.write('\n')
	file.close()



	##======== timing
	time_end = time.time()
	print "time spent on this gene is",
	print time_end - time_start



