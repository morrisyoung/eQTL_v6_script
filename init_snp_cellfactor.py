## tasks:
##	1. prepare the dosage matrix, transpose, inverse;
##	2. load the factor matrix, and multiply by the above matrix, to get the initialization of the snp_cellfactor parameter
## NOTE:
## issues!!!
##	1. (May.25) can't load all the SNPs for all samples, not to say do the regression (96G memory)
## solution!!!
##	1. do the single-point regression for all the SNPs; but this will lose some dependency among the SNPs --> see #2
##	2. split the SNP space into several sub-space, and do regression for each; we lose some dependency, but we have no other solutions beyond this

## This script is for step#2




import numpy as np



num_snp = 0





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





if __name__ == "__main__":


	##=====================================================
	##============ prepare the dosage matrix ==============
	##=====================================================
	# under the below folder we will have the dosage data grouped by chromosomes
	filename = "../../GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/genotype_imputed/genotype_450_dosage_matrix_qc"
	# we need two:
	# "/chrX/SNP_dosage_ID.txt" and "/chrX/SNP_info.txt"

	#==== get the total number of SNPs first
	num_snp = 0
	for i in range(1, 23):
		chr = str(i)
		file = open(filename + "/chr" + chr + "/SNP_info.txt", 'r')

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			num_snp += 1
		file.close()

	print "there are # total number of SNPs:",
	print num_snp







	##=====================================================
	##============ solving the linear sysmte ==============
	##=====================================================
	# reference:
	#	http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq
	#	http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html#solving-linear-least-squares-problems-and-pseudo-inverses

	## TODO:
	part = 1


	#==== dosage
	Data = np.load("../data_processed/init_dosage_train_m" + str(part) + ".npy")
	print len(Data)
	print len(Data[0])

	#==== factor (and factor_coef)
	factor = np.load("../data_processed/nmf_factor.npy")
	factor_coef = np.load("../data_processed/nmf_coefficient.npy")
	# tune down factor power, and tune up factor_coef power
	max = np.amax(factor)
	min = np.amin(factor)

	#scaler = xxx	# NOTE: this doesn't work!!!

	factor = factor * scaler
	factor_coef = factor_coef * (1/scaler)

	np.save("../data_processed/nmf_factor_rescale", factor)
	np.save("../data_processed/nmf_coefficient_rescale", factor_coef)


	#reverse logistic








	## won't consider concatenate, as the pinv need larger than 96G memory
	'''
	Data = np.concatenate((Data, Data), axis=1)
	print len(Data)
	print len(Data[0])

	Data_inv = np.linalg.pinv(Data)
	print len(Data_inv)
	print len(Data_inv[0])
	'''


	## TODO: I will do both lstsq and pinv --> it seems both two work, but the second takes less memory (job IDs: 483953; 483955)
	##==== inv
	'''
	Data_inv = np.linalg.pinv(Data)
	print len(Data_inv)
	print len(Data_inv[0])

	factor = np.random.rand(5170, 400)

	m = np.dot(Data_inv, factor)

	print len(m)
	print len(m[0])
	'''


	##==== lstsq
	#factor = np.random.rand(5170, 400)
	m = np.linalg.lstsq(Data, factor)[0]
	print len(m)
	print len(m[0])

	np.save("../data_processed/init_dosage_para_m" + str(part), m)



