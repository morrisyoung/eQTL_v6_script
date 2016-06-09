## given the Pearson correlation results (between prediction and real data), plot them for all genes

## input of this program:
##	1. correlation values ("../analy_result/analy_corr_tissueX.txt")
##	2. gene list in order ("../analy_result/analy_gene_list_tissueX.npy")




import numpy as np
import matplotlib.pyplot as plt






# information table (colors)
color_map = {'1':'m', '2':'#81b1d2', '3':'#ffed6f', '4':'r', '5':'#EEEEEE', '6':'#cbcbcb', '7':'#6d904f', '8':'y', '9':'#E24A33', '10':'#0072B2', '11':'#f0f0f0', '12':'0.40', '13':'blue', '14':'#fc4f30', '15':'#bfbbd9', '16':'#ccebc4', '17':'c', '18':'#A60628', '19':'#988ED5', '20':'g', '21':'#bcbcbc', '22':'#FFB5B8', 'X':'m', 'Y':'r', 'MT':'blue'}






if __name__ == "__main__":

	###
	tissue_rep = {}
	file = open("../analy_result/tissue_index_map.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		index = int(line[1])
		tissue_rep[tissue] = index
	file.close()

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



	########
	######## tissue loop
	########
	for tissue in tissue_rep:

		##===================================================== corr_rep and gene_list =====================================================
		corr_rep = {}
		filename = "../analy_result/analy_corr_tissue" + str(tissue_rep[tissue]) + ".txt"
		file = open(filename, 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]
			corr = float(line[1])
			corr_rep[gene] = corr
		file.close()
		
		filename = "../analy_result/analy_gene_list_tissue" + str(tissue_rep[tissue])
		gene_list = np.load()



		##===================================================== plot =====================================================
		# with gene_list and corr_rep
		plt.figure(1)
		for gene in gene_list:
			corr = corr_rep[gene]

			chr = int(gene_tss[gene][0])
			color = color_table[chr-1]
			plt.plot(i, corr, color, marker = 'o', alpha=0.7)

			'''
			## add the gene id if the corr is high enough
			if corr >= 0.5:
				print gene,
				print corr
			'''

		#plt.axis([0, 20000, -1, 1])
		plt.xlabel('Expressed genes (coding and non-coding) from all 22 chromosomes')
		plt.ylabel('Pearson correlation of gene expression level')
		plt.title('Model testing for the multi-linear regression of cis- SNPs (+-1Mb)')
		plt.title('Model testing (' + tissue + ') for the multi-linear regression of cis- SNPs (+-1Mb)')
		plt.grid(True)


		plt.show()




