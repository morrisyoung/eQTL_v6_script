## test the prediction power of the cis- SNPs (multivariate linear regression model)

## TODO: I will try different data pre-processing methods


## notes:
##	1. when multiple tissues are involved, it makes no sense to predict the gene expression just by the cis- regulators; we must have tissue indicators as a predictor in this setting, and do in a tissue-specific fashion --> this is just equal to do linear regression for each tissue separately (the intercept variables will expalin tissue effects!!!)  otherwise the same predictor (the genotype of the same individual) will have multiple predictions for different tissues; although there are only 75+ training samples to do this regression, we need to do that
##	2. from this experiment, we want to see that the simple combination of cis- regulators and tissue effects don't have much predictive performance; and that's why we try to build a better model --> two directions: a. hierarchical linear regression model; b. our neural net model
##	3. [TODO] for the neural model, how can we let parameters share information across different tissues? we probably need to add a prior on top of the model





import time
import numpy as np
import scipy.linalg	# NOTE: the way to import (specially for Scipy)



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







## NOTE: I'm testing various normalization for the data








##==================
##==== sub-routines
##==================
# get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
def sample_to_individual(s):
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
				snp_dosage_rep[individual][-1].append(dosage)
			file.close()
	return




if __name__ == "__main__":

	print "working on the multi-linear regression of cis- SNPs for all genes"



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
	print "there are # of individuals:",
	print num_individual
	snp_dosage_load()


	#DEBUG
	for individual in snp_dosage_rep:
		print "total # of SNPs:",
		count = 0
		for i in range(22):
			count += len(snp_dosage_rep[individual][i])
		print count
		break



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
			snp_pos_list[-1].append(pos)
		file.close()





	########
	######## tissue loops
	########
	tissue_rep = {}
	file = open("../data_processed/list_samples_train.txt", 'r')
	tissue_count = 1
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		tissue_rep[tissue] = tissue_count
		tissue_count += 1
	file.close()
	file = open("../analy_result/tissue_index_map.txt", 'w')
	for tissue in tissue_rep:
		file.write(tissue + '\t' + str(tissue_rep[tissue]) + '\n')
	file.close()

	for tissue in tissue_rep:


		print "working on tissue:",
		print tissue
		time_start = time.time()





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
			if line[0] == tissue:
				sample_list = line[1:]
				for i in range(len(sample_list)):
					sample = sample_list[i]
					sample_rep[sample] = []

		file.close()


		##=================== query the RPKM according to the above list
		print "get expression matrix for these training samples..."
		#file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')
		file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_z", 'r')
		#file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_Gaussian_rank", 'r')

		###
		sample_list = ((file.readline()).strip()).split('\t')[1:]
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
			rpkm_list = map(lambda x: float(x), line[1:])

			for i in range(len(rpkm_list)):
				rpkm = rpkm_list[i]
				if i in index_rep:
					sample = index_rep[i]
					sample_rep[sample].append(rpkm)
		file.close()

		print "there are # of samples for training:",
		print len(sample_rep)

		###
		gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
		for i in range(len(gene_list)):
			gene = gene_list[i]
			gene_index_map[gene] = i

		###
		num_gene = len(gene_list)
		print "number of genes:",
		print num_gene

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
		'''
		# gene_cis_index	NOTE: only need to do the processing routine once (for the same gene annotation, and the same set of SNPs)
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
				if flag1 == 1 and flag2 == 1:
					break
			gene_cis_index[gene] = (start, end)
		file = open("../data_processed/gene_cis_range.txt", 'w')
		for gene in gene_cis_index:
			start = gene_cis_index[gene][0]
			end = gene_cis_index[gene][1]
			file.write(gene + '\t' + str(start) + '\t' + str(end) + '\n')
		file.close()
		'''


		##==== for future I won't need to do above again
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
				individual = sample_to_individual(sample)
				genotype_matrix.append([])
				for k in range(start, end+1):
					dosage = snp_dosage_rep[individual][chr-1][k]
					genotype_matrix[-1].append(dosage)
				genotype_matrix[-1].append(1)  # NOTE: we need the intercept, and we need to manually add the input item (1) for the intercept
			genotype_matrix = np.array(genotype_matrix)
			## sample:
			#X = np.array([[1,2,3,1], [2,4,6,1], [3,6,9,1]])  # 1x1 + 2x2 + 3x3 + 1
			#y = np.array([14, 28, 42])
			#y = np.array([15, 29, 44])
			#m = np.linalg.lstsq(X, y)[0]
			try:
				## NOTE: the following gene can't make SVD converge:
				##	"ENSG00000109790.12"
				m = np.linalg.lstsq(genotype_matrix, expression_array)[0]
				## NOTE: try Scipy:
				#m = scipy.linalg.lstsq(genotype_matrix, expression_array)[0]
				para_rep[gene] = np.array(m)  ## there is an extra intercept here!!!
			#except ValueError:
			except:
				print "error for gene:",
				print gene
				print "genotype_matrix is:"
				print genotype_matrix
				print "expression_array is:"
				print expression_array
				## write the matrix into a file
				np.save("./temp/" + gene + "_geno", genotype_matrix)
				np.save("./temp/" + gene + "_gene", expression_array)
				'''
				file = open("./temp/" + gene + ".txt", 'w')
				for i in range(len(genotype_matrix)):
					for j in range(len(genotype_matrix[i])):
						dosage = genotype_matrix[i][j]
						file.write(str(dosage) + '\t')
					file.write('\n')
				file.close()
				'''

				## NOTE: in this case, I will instead inverse the matrix, and then dot product to get the coefficients
				## NOTE: for some, still SVD can't succeed --> SO: if can't do, don't do
				'''
				genotype_matrix_inv = np.linalg.pinv(genotype_matrix)
				m = np.dot(genotype_matrix_inv, expression_array)
				para_rep[gene] = m
				'''



		##===================================================== save all learned parameters ===================================================
		filename = "../analy_result/analy_para_tissue" + str(tissue_rep[tissue]) + ".txt"
		file = open(filename, 'w')
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
		print "time spent on training this tissue is",
		print time_end - time_start








		##======================================================================================================================
		##======================================================================================================================
		## testing part (for each tissue)
		##======================================================================================================================
		##======================================================================================================================








		print "testing the model, output the Pearson correlation (predicted values and real values)"
		time_start = time.time()


		##===================================================== expression (testing) =====================================================
		##================== get the current training tissues and samples inside (list)
		print "get testing samples..."
		file = open("../data_processed/list_samples_test.txt", 'r')
		sample_rep = {}
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			if line[0] == tissue:		## only need this tissue
				sample_list = line[1:]
				for i in range(len(sample_list)):
					sample = sample_list[i]
					sample_rep[sample] = []

		file.close()


		##=================== query the RPKM according to the above list
		print "get expression matrix for these testing samples..."
		#file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')
		file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_z", 'r')
		#file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_Gaussian_rank", 'r')

		###
		sample_list = ((file.readline()).strip()).split('\t')[1:]
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
			rpkm_list = map(lambda x: float(x), line[1:])

			for i in range(len(rpkm_list)):
				rpkm = rpkm_list[i]
				if i in index_rep:
					sample = index_rep[i]
					sample_rep[sample].append(rpkm)
		file.close()

		print "there are # of samples for testing:",
		print len(sample_rep)

		###
		gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
		for i in range(len(gene_list)):
			gene = gene_list[i]
			gene_index_map[gene] = i

		###
		num_gene = len(gene_list)
		print "number of genes:",
		print num_gene

		### this is crucial for building the expression matrix, as we need the order of all samples
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
			chr = line[1]  # don't int() here, as there are XYMT genes
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


		##================================================== testing all genes, get the corr_rep ==============================================
		corr_rep = {}
		for i in range(len(gene_list)):

			gene = gene_list[i]

			if gene not in para_rep:
				continue

			chr = int(gene_tss[gene][0])
			start = gene_cis_index[gene][0]
			end = gene_cis_index[gene][1]

			##
			expression_array_real = []
			for j in range(len(sample_list)):
				rpkm = expression_matrix[j][i]
				expression_array_real.append(rpkm)
			expression_array_real = np.array(expression_array_real)


			##
			expression_array_exp = []
			for j in range(len(sample_list)):
				genotype_array = []
				sample = sample_list[j]
				individual = sample_to_individual(sample)
				for k in range(start, end+1):
					dosage = snp_dosage_rep[individual][chr-1][k]
					genotype_array.append(dosage)
				genotype_array.append(1)  # we need the intercept
				genotype_array = np.array(genotype_array)
				## expected expression level
				rpkm = np.inner(genotype_array, para_rep[gene])
				expression_array_exp.append(rpkm)
			expression_array_exp = np.array(expression_array_exp)



			## DEBUG
			filename = "../analy_z_result/gene_" + gene
			np.save(filename, expression_array_exp)



			##
			corr = np.corrcoef(expression_array_real, expression_array_exp)[0][1]
			corr_rep[gene] = corr




		##===================================================== save corr_rep =====================================================
		filename = "../analy_result/analy_corr_tissue" + str(tissue_rep[tissue]) + ".txt"
		file = open(filename, 'w')
		for gene in corr_rep:
			file.write(gene + '\t' + str(corr_rep[gene]) + '\n')
		file.close()

		gene_list_new = []
		# save gene_list in order:
		for gene in gene_list:
			if gene in para_rep:
				gene_list_new.append(gene)
		gene_list_new = np.array(gene_list_new)
		filename = "../analy_result/analy_gene_list_tissue" + str(tissue_rep[tissue])
		np.save(filename, gene_list_new)




		##======== timing
		time_end = time.time()
		print "time spent on this tissue is",
		print time_end - time_start






		## DEBUG:
		break








	print "done!..."



