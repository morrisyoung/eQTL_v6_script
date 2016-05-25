## initialize the cellenv variables (and coefficients) with PCA (if necessary, probably also with Factor Analysis method)
## data dependency:
##	../data_processed/list_samples_train.txt
##	../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize
## several steps:
##	get the training samples
##	perform PCA
##	visualize




import numpy as np
from sklearn.decomposition import PCA
#import matplotlib.pyplot as plt
from sklearn.decomposition import NMF




n_factor = 400		# TODO




if __name__ == "__main__":


	##============================================
	##============= processing data ==============
	##============================================
	'''
	print "cellenv factor init..."

	##==== load all the training samples (IDs)
	sample_rep = {}
	file = open("../data_processed/list_samples_train.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		for sample in line[1:]:
			sample_rep[sample] = 1
	file.close()




	##==== load and process the expression matrix
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')

	# get the sample list and index rep
	sample_list = []
	index_rep = {}
	line = (file.readline()).strip()
	line = (line.split('\t'))[1:]	# NOTE here, we pruned the first one
	for i in range(len(line)):
		sample = line[i]
		if sample in sample_rep:
			sample_list.append(sample)
			index_rep[i] = 1
	Data = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		Data.append([])

		line = (line.split('\t'))[1:]
		for i in range(len(line)):
			if i in index_rep:
				value = float(line[i])
				Data[-1].append(value)
	Data = np.array(Data)
	file.close()
	print "loaded data..."


	##==== saving the processed expression matrix
	np.save("../data_processed/expression_train", Data)
	np.save("../data_processed/expression_train_samples", np.array(sample_list))
	print "saving done"

	'''



	##================================
	##============= PCA ==============
	##================================
	'''
	##==== saving the processed expression matrix
	print "loading data..."
	Data = np.load("../data_processed/expression_train.npy")
	sample_list = np.load("../data_processed/expression_train_samples.npy")



	##==== PCA
	## the matrix is (gene x sample)
	print "performing PCA..."
	pca = PCA(n_components=n_factor)
	pca.fit(Data)
	Y2 = pca.components_
	Y1 = pca.transform(Data)
	variance = pca.explained_variance_ratio_

	print variance


	# DEBUG
	print len(Y1),
	print len(Y1[0])
	print len(Y2),
	print len(Y2[0])



	##==== save PCA results (two matrices, for coefficient matrix and factor matrix; and also the sample_list)
	np.save("../data_processed/pca_coefficient", Y1)
	np.save("../data_processed/pca_factor", Y2)
	np.save("../data_processed/pca_variance", variance)
	print "saving done"
	
	'''





	##================================
	##============= NMF ==============
	##================================
	## link:
	##	http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html#sklearn.decomposition.NMF.fit_transform
	##==== saving the processed expression matrix
	print "loading data..."
	Data = np.load("../data_processed/expression_train.npy")
	sample_list = np.load("../data_processed/expression_train_samples.npy")


	##==== NMF
	print "performing NMF..."
	model = NMF(n_components=n_factor, init='random', random_state=0, max_iter=1000)
	Y1 = model.fit_transform(Data)
	Y2 = model.components_
	iteration = model.n_iter_
	error = model.reconstruction_err_


	# DEBUG
	print iteration
	print error
	print len(Y1),
	print len(Y1[0])
	print len(Y2),
	print len(Y2[0])



	##==== save NMF results (two matrices, for coefficient matrix and factor matrix)
	np.save("../data_processed/nmf_coefficient", Y1)
	np.save("../data_processed/nmf_factor", Y2)
	print "saving done"


	'''
	## visualized saving
	Y1 = np.load("../data_processed/nmf_coefficient.npy")
	Y2 = np.load("../data_processed/nmf_factor.npy")
	file = open("../data_processed/nmf_coefficient.txt", 'w')
	for i in range(len(Y1)):
		for j in range(len(Y1[i])):
			value = str(Y1[i][j])
			file.write(value + '\t')
		file.write('\n')
	file.close()

	file = open("../data_processed/nmf_factor.txt", 'w')
	for i in range(len(Y2)):
		for j in range(len(Y2[i])):
			value = str(Y2[i][j])
			file.write(value + '\t')
		file.write('\n')
	file.close()
	'''






	
	'''
	##======================================
	##============= visualize ==============
	##======================================
	Y1 = np.load("../data_processed/pca_coefficient.npy")
	Y2 = np.load("../data_processed/pca_factor.npy")
	variance = np.load("../data_processed/pca_variance.npy")



	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	#image = np.random.poisson(10., (100, 80))
	#i = ax.imshow(image, interpolation='nearest')
	i = ax.imshow(Y1[:500], interpolation='nearest')
	fig.colorbar(i)  # note that colorbar is a method of the figure, not the axes
	plt.show()



	## visualized saving
	Y1 = Y1 / 10
	file = open("../data_processed/pca_coefficient.txt", 'w')
	for i in range(len(Y1)):
		for j in range(len(Y1[i])):
			value = str(Y1[i][j])
			file.write(value + '\t')
		file.write('\n')
	file.close()

	Y2 = Y2 * 10
	file = open("../data_processed/pca_factor.txt", 'w')
	for i in range(len(Y2)):
		for j in range(len(Y2[i])):
			value = str(Y2[i][j])
			file.write(value + '\t')
		file.write('\n')
	file.close()

	file = open("../data_processed/pca_variance.txt", 'w')
	for i in range(len(variance)):
		value = str(variance[i])
		file.write(value + '\n')
	file.close()
	'''




