# RNAseq_vc
## How to call variants on references using RNA-seq data

You just need a reference on which to call the variants and libraries for each sample, and a *list.sample* file is recommended to speed up some below commands.
In case of many sample per species, you would expect to have the same variants(intra-population) or not!

You can both generate a consensus alternative reference of your species, or have as many alternative reference as the number of sample you have.

## 1. MAPPING RNA-seq LIBRARIES ON REFERENCE
Use command in *bcftools_vc.sh* to generate $sample'.alternative.calls.tab' files, in which you have the variant calling format only for alternative sites.
Change the script accordingly with your needs.

## 2. COMPARE ALTERNATIVE SITES AMONG SAMPLES
In this way you can decide to generate a file of variant calling sites shared by all of you sample, if your aim is to have a consensus sequence, and you know that the variance across your sample is low (intra-population) and the alternative calls that is not shared among samples could be 1. a wrong call in a region of low coverage, or 2. a rare variant (in case of many samples and just one or a few having the alternative site). So the threshold will be equal to the number of samples you have.
```
threshold=$(wc -l list.sample | awk '{print $1}')
python gpt3_compare.py $threshold $PWD > 'out_of_'$threshold'.output.vcalling.txt'
```
Or maybe you want to have a file in which for each sites that differes from the reference you know what base has been called on that sites. To know if there is a pattern among your samples (maybe sub-populations?). Store all the '.alternative.calls.tab' files in a directory and just run...
```
python gpt4_compare.py $all_alternative_dir
sed -i 's/,/\t/g' merged_alternative_values.csv
```
Moreover, you may want to call each variant indipendently and than use those files to fix the reference for each sample, in order to have many alternative sequences. Just do something like this...
```
for sample in $(cat list.sample)
  do mkdir $sample
  mv $sample'.alternative.calls.tab' $sample/.
  python gpt3_compare.py $sample $PWD/$sample > 'out_of_'$sample'.output.vcalling.txt'
  done
```
## 3. CORRECT THE REFERENCE
One more steps to go.Now we need to fix the refernce using the alternative calls previously generated. If look for the unique consensus sequence, just use the python command with the number of threshold as $sample var.
The reference fasta must be the same used at the beginning to map the RNA-seq libraries (just a matter of coordinates).
```
for sample in $(cat list.sample)
  do python correct.py $sample $PWD $reference_fasta
  done
```
## EXPRESSION PLOTTING
You may be interested in plotting the expression of your RNA-seq data along your reference, so to spot immediately region of low coverage.
Use the following!
```
samtools depth $sample'.sorted.bam' > $sample'.coverage.txt'
python plot.coverage.py
```
