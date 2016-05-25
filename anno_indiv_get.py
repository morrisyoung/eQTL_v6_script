## this is used to get the individual list from some source file.




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

	file = open("../data_source/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_sample_IDs.txt", 'r')
	file1 = open("../data_processed/list_individuals.txt", 'w')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		sample = line
		individual = sample_to_individual(sample)
		file1.write(individual + '\n')

	file.close()
	file1.close()

