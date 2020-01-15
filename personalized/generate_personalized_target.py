import subprocess


def read_target_region(filename):
    """
    reads a fasta file and generates the target region as a string
    :param filename: fasta file (target region)
    :return: string without the line containing '>'
    """
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content][1:]
    return ''.join(content)


def detect_variant(target_filename, reads_filename):
    """
	Reads the target, aligns the seq reads on the target using bowtie2, uses Pilon to detect the genetic variations, then returns the improved target sequence.
    :param target_filename: the target region filename as a fasta
    :param reads_filename: the fastq/fasta reads file. Currently, only unpaired
    :return: returns the target region after detecting personalized variations as a string
	"""
    with open('personalized/commands_py.sh', 'w') as f:
        f.write('cd ./personalized\n')
        f.write('bowtie2-build -f ../' + target_filename + ' test\n')
        f.write('bowtie2 --local -x test -U ../' + reads_filename + ' -S sam_out.sam\n')
        f.write('samtools view -bS sam_out.sam > bam_out.bam\n')
        f.write('samtools sort bam_out.bam > bam_sorted.bam\n')
        f.write('samtools index bam_sorted.bam\n')
        f.write('java -Xmx16G -jar pilon-1.23.jar --genome ../' + target_filename + ' --unpaired bam_sorted.bam\n')
        f.close()
    return_code = subprocess.call(['bash', 'personalized/commands_py.sh'])
    if return_code != 0:
        print('Something went wrong!')
        exit(-1)
    variant_target = read_target_region('personalized/pilon.fasta')
    return variant_target


if __name__ == '__main__':
    variant = detect_variant('inputs/target20k.txt', 'inputs/read_staphylo.fastq')
    original_target = read_target_region('inputs/target20k.txt')
    print(len(variant))
    print(len(original_target))