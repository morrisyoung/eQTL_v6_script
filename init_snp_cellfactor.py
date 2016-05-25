## tasks:
##	1. prepare the dosage matrix, transpose, inverse;
##	2. load the factor matrix, and multiply by the above matrix, to get the initialization of the snp_cellfactor parameter


import numpy as np




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


	## NOTE: (May.25) can't load (96G memory)


	# under the below folder we will have the dosage data grouped by chromosomes
	filename = "../../GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/genotype_imputed/genotype_450_dosage_matrix_qc"
	# we need two:
	# "/chrX/SNP_dosage_ID.txt" and "/chrX/SNP_info.txt"


	Data = []
	sample_list = np.load("../data_processed/expression_train_samples.npy")

	count = 0

	for sample in sample_list:
		print count,
		count += 1
		print sample


		individual = sample_to_individual(sample)

		Data.append([])
		for i in range(1, 23):
			chr = str(i)
			file = open(filename + "/chr" + chr + "/SNP_dosage_" + individual + ".txt", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				dosage = float(line)
				Data[-1].append(dosage)
			file.close()
	print "now loaded"

	Data = np.array(Data)

	print "now nped"

	np.transpose(Data)

	print "now transposed"



	print len(Data)
	print len(Data[0])

	np.save("../data_processed/init_dosage_train", Data)






	##=====================================================
	##============ solving the linear sysmte ==============
	##=====================================================
	# reference:
	#	http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq
	#	http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html#solving-linear-least-squares-problems-and-pseudo-inverses



