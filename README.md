# eQTL_v6_script

I use these scripts to process GTEx.v.6 data.


## 1. Gene expression data

1. sample\_tissue\_preprocess.py

## 2. Genotype data

In the imputed data directory ("phg000520.v2.GTEx\_MidPoint\_Imputation.genotype-calls-vcf.c1"), we have the VCF-format genotype data. Create the following folder:

1. ./genotype\_imputed
2. ./genotype\_imputed/genotype\_vcf
3. ./genotype\_imputed/genotype\_processed

Then run the two sections (separately) of "genotype\_vcf\_parser.py" under "./genotype\_imputed/" and for each chromosome separately (as the disk quota issue), to perform the following two functions:

1. splitting the VCF file into header file and sub-VCF file for each chromosome (saved to "./genotype\_imputed/genotype\_vcf/")
2. extracting the tped/dosage/snp\_exclusion\_list information from each sub-VCF file for each chromosome (saved to "./genotype\_imputed/genotype\_processed/")

We can find the tfam file under "../phg000520.v2.GTEx\_MidPoint.genotype-qc.MULTI/Imputation/".

After that, we can use the procedure similar to GTEx.v.4 to process the data (snp QC exclusion, pruning; here is the [link](https://github.com/morrisyoung/eQTL_v4_script#5-the-pipeline-for-genotype-qc-and-ld-pruning)).

**Note:** We need "chrX.tped" and "chrX.tfam" (--tfile) to perform the snp QC exclusion (witn snp\_exclusion\_list), and that will generate bed/fam/bim files that we can do the LD pruning. And we can finally extract the dosage information from "chrX.dosage".
