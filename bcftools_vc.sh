for sample in $(cat list.sample)
do

cd $sample

reference_fasta=BaAt_mit.fasta
reference=${reference_fasta%.fasta}

#bowtie2-build -f $reference_fasta $reference
bowtie2 --threads 12 -x ../$reference -q -1 *_1.trimmed.fp.fastq.gz -2 *_2.trimmed.fp.fastq.gz -S $sample'.sam'
#samtools view --threads 12 -bS $1'.sam' > $1'.bam'
samtools sort --threads 12 $sample'.sam' -o $sample'.sorted.bam'

bcftools mpileup --threads 12 --max-depth 5000 -Ou -f ../$reference_fasta $sample'.sorted.bam' | bcftools call -c --pval-threshold 0.25 -Oz -o calls.vcf.gz
gzip -d calls.vcf.gz

#set the spacing for the next cat
IFS='
'

for line in $(cat calls.vcf | grep -v '^#')
    do value=$(echo $line | awk -F '\t' '{print $5}')
    if [ $value != "." ]
        then echo $line >> $sample'.alternative.calls.tab'
        fi
    done

cd ..
