## to partition the individuals (NOTE: other than samples) into training set and testing set (75% and 25%), and organize the files by tissues (as required by the training program).
## input file:
##	phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_100 #TODO: the threshold


import numpy as np



filter = 100	# TODO
ratio = 0.75	# TODO



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




	##==== load all samples
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + str(filter), 'r')

	tissue_sample_rep = {}		# element as list of samples (from that tissue)
	sample_tissue_map = {}		# element as tissue type (of that sample)
	individual_sample_rep = {}	# element as list of samples (from that individual)

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]

		## tissue_sample_rep
		tissue_sample_rep[tissue] = sample_list

		## sample_tissue_map
		for i in range(len(sample_list)):
			sample = sample_list[i]
			sample_tissue_map[sample] = tissue

		## individual_sample_rep
		for i in range(len(sample_list)):
			sample = sample_list[i]
			individual = sample_to_individual(sample)
			if individual not in individual_sample_rep:
				individual_sample_rep[individual] = [sample]
			else:
				individual_sample_rep[individual].append(sample)

	file.close()




	##==== split
	## targets:
	##	list_samples_train.txt
	##	list_samples_test.txt
	## pick up the training individuals and testing individuals, and save them
	file1 = open("../data_processed/list_samples_train.txt", 'w')
	file2 = open("../data_processed/list_samples_test.txt", 'w')

	individual_list = []
	for individual in individual_sample_rep:
		individual_list.append(individual)
	arr = np.arange(len(individual_list))
	np.random.shuffle(arr)


	## the following script is for saving the samples indexing by individual IDs
	'''
	i = 0
	while i < len(individual_list) * ratio:
		index = arr[i]
		individual = individual_list[index]
		file1.write(individual + '\t')
		## write the samples from this individual, to training set
		for j in range(len(individual_sample_rep[individual])):
			sample = individual_sample_rep[individual][j]
			file1.write(sample + '\t')
		file1.write('\n')
		i += 1


	while i < len(individual_list):
		index = arr[i]
		individual = individual_list[index]
		file2.write(individual + '\t')
		## write the samples from this individual, to testing set
		for j in range(len(individual_sample_rep[individual])):
			sample = individual_sample_rep[individual][j]
			file2.write(sample + '\t')
		file2.write('\n')
		i += 1
	'''


	## we should save indexing by tissues
	rep_train = {}
	rep_test = {}
	for tissue in tissue_sample_rep:
		rep_train[tissue] = []
		rep_test[tissue] = []

	## splitting the individuals
	i = 0
	while i < len(individual_list) * ratio:
		index = arr[i]
		individual = individual_list[index]
		for sample in individual_sample_rep[individual]:
			tissue = sample_tissue_map[sample]
			rep_train[tissue].append(sample)
		i += 1
	while i < len(individual_list):
		index = arr[i]
		individual = individual_list[index]
		for sample in individual_sample_rep[individual]:
			tissue = sample_tissue_map[sample]
			rep_test[tissue].append(sample)
		i += 1

	## saving the samples by tissues
	for tissue in rep_train:
		file1.write(tissue + '\t')
		for sample in rep_train[tissue]:
			file1.write(sample + '\t')
		file1.write('\n')
	for tissue in rep_test:
		file2.write(tissue + '\t')
		for sample in rep_test[tissue]:
			file2.write(sample + '\t')
		file2.write('\n')

	file1.close()
	file2.close()




