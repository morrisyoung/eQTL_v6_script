## tasks:
##	do the linear regression for the cis- SNPs and the gene expression






if __name__ == "__main__":



	##=============================
	##========== tss ==============
	##=============================
	## get the chr and tss for all genes (including X, Y and MT genes)
	file = open("../data_source/gencode.v19.genes.patched_contigs.gtf", 'r')
	file1 = open("../data_processed/gene_tss.txt", 'w')

	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split()
		chr = line[0]
		tss = line[3]
		gene = line[9][1: -2]

		type = line[2]
		if type == 'transcript':
			file1.write(gene + '\t' + chr + '\t' + tss + '\n')

	file.close()
	file1.close()



	##=======================================
	##========== X, Y, MT list ==============
	##=======================================
	## get the list of all X, Y, MT genes
	file = open("../data_processed/gene_tss.txt", 'r')
	file1 = open("../data_processed/gene_xymt.txt", 'w')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		tss = line[2]

		if chr == 'X' or chr == 'Y' or chr == 'MT':
			file1.write(gene + '\n')
	file.close()
	file1.close()



