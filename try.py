## this is processing script





if __name__ == "__main__":


	print "hello world."




	## check the max and min in the transformed expression profile
	"""
	max = 1
	min = 1

	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')
	file.readline()
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		rpkm_list = map(lambda x: float(x), line[1:])
		for rpkm in rpkm_list:
			if rpkm > max:
				max = rpkm
			if rpkm < min:
				min = rpkm


	file.close()

	print max
	print min
	"""
