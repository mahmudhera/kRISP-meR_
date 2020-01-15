bowtie2-build -f target20k.txt test
bowtie2 --local -x test -U read_staphylo.fastq -S sam_out.sam
samtools view -bS sam_out.sam > bam_out.bam
samtools sort bam_out.bam > bam_sorted.bam
samtools index bam_sorted.bam
java -Xmx16G -jar pilon-1.23.jar --genome target20k.txt --unpaired bam_sorted.bam

def generate():
    open('commands.sh')