![image](https://img.shields.io/badge/%20-linux-orange)
![image](https://img.shields.io/badge/%20-python-blue)
![image](https://img.shields.io/badge/crispr-referencefree-yellowgreen)
# kRISP-mER
Reference free guide RNA designing tool for CRISPR

## What kRISP-mER is
This a tool to generate personalized guide RNAs for CRISPR without using a reference genome. Instead, kRISP-mER works using the sequenced reads, and a genomic target location (location where CRISPR cleavage is intended).

This tool is designed for:
1. Linux
1. Python
1. If on-target-activity scores are required, then particularly, **Python 2.7**

Other dependencies are elaborated separately.

## Dependencies you need to install
The following installation instructions are _only to help you out_. These installation instructions are **NOT** mandatory to follow. You can install these any way you like. However, if you are having a hard time doing so by yourself, you may find these instructions useful.
* **samtools**: You can install samtools using the following commands:
```shell script
sudo apt-get update -y
sudo apt-get install -y samtools
```
* **Biopython**: install using: `pip install biopython`
* **Python binding of Jellyfish** (gmarics project). This is quite tricky. Need to assess this in detail later. Tried to do the following:
```shell script
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
```
Python binding means that if you try to import jellyfish from a python script (the code is: `import jellyfish`), that will work and you will be able to invoke Jellyfish program from python.
If the installation steps mentioned above does not bind with python, (although the documentation of Jellyfish does say that this should): you may try the swig binding instructions from https://github.com/gmarcais/Jellyfish/blob/master/swig/Readme.md).
```shell script
cd swig/python
python setup.py build
python setup.py install --prefix=$PREFIX
```
* **Bowtie2**: You can install this with With `Bioconda`. With `Bioconda` installed, you should be able to install Bowtie 2 with `conda install bowtie2`. Details are found here: `http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2`
* **Java**: kRISP-mER uses Pilon to find a better sequence after seq reads are aligned to a base sequence. Pilon is run as a plain jar file. In order to check this requirement, you may want to test `java --version` and `jre --version`
* **numpy**: Simple installation with pip: `pip install numpy`. You may also work with some package manager, such as anaconda, in which case numpy should come built in.
* **scipy**: Simple installation with pip: `pip install scipy`. You may also work with some package manager, such as anaconda, in which case scipy should come built in.
* **sklearn**: If the package manager you are using does not already have scikit-learn installed, you can install using `pip install scikit-learn==0.16.1` (This very specific version is important to determine on target activity scores)
* **pickle**: This python package is required to determine on-target-activity as well. This should already be installed in python 2 and 3. If not, you need to manually install this.

## How to run
Once you have the dependencies installed, running kRISP-mER is easy. You need to:
1. **Download** the github repository
1. **Locate** your sequenced reads file (can be FASTA or FASTQ) and the target-region file (_must_ be FASTA); these files can be anywhere on the filesystem
1. **Run** the tool with the command: `python krispmer.py <reads_filename> <target_filename> <num_mismatches> <output_file>` _after navigating to the directory containing the file krispmer.py_
The arguments, along with other options can be seen using `python krispmer.py -h`.

**Note**: You may be limited by your workstation resources. kRISP-mER works on sequenced reads by counting the k-mers using Jellyfish. If the reads file is very large, Jellyfish may be out of memory and kill the program. If that happens, you may want to increase the resources you have. kRISP-mER stores the Jellyfish binary file and uses that to calculate the gRNA scores. Therefore, you have to make sure that the workstation filesystem has enough free storage to store this binary file (whose size depends on the size of the sequenced reads). 

## Options available in kRISP-mER
Usage:
```shell script
python krispmer.py [-h] [-p] [-e] [-s] [-v] [-t] [-n] [-c CUTOFF_SCORE] [-g GENOME] [-a ALTPAMS [ALTPAMS ...]]
                   <reads_file> <target_file> <output_file> <max_hd>
```
kRISP-mER allows you to design guide RNAs with WGS shotgun reads (in a FASTA or FASTQ file), and a target-region (a FASTA file). With these two, you also have to tell the program the number of mismatches to consider when designing a gRNA. kRISP-mER allows upto 3 mismatches. kRISP-mER does not consider indels (like other established gRNA designing tools). You also have to tell the program the name of the output csv file, where the gRNAs along with their inverted specificity scores and strand information is to be stored.

Besides these four positional (mandatory) arguments, you can also do the following.
1. `-h`: You can see help with `-h` flag
1. `-p`: You can tell the program to let you manually choose the read coverage with the flag `-p`. By default, the program automatically calculates the read coverage from the k-spectrum of the sequenced reads
1. `-e`: You can specify the program to calculate the prior probabilities with EM algorithm on the k-spectrum of the sequenced reads with `-e` flag. See the full paper to understand what the EM algorithm does. By default, EM does not work.
1. `-s`: You can choose to exclude the gRNAs that contain a stop codon with the flag `-s`. By default, these gRNAs are included in the final output list
1. `-v`: You can tell kRISP-mER to detect the genetic variation of the individual in the target site using the flag `-v` (which is done with a combination of bowtie2, samtools and Pilon). By default, the target is used as is.
1. `-n`: You can specify to detect guides from only the positive (5'-3') strand with the flag `-n`. By default, both strands are considered.
1.  `-c INT`: You can pass a cut-off score for the guides with the flag `-c`. kRISP-mER will then drop all the guides with inverted-specificity higher than the cut-off.
1. `-a PAM1 PAM2 ...`: You can provide kRISP-mER with a list of PAMs to consider with `-a` flag. By default, NGG PAMs are considered.

## Do not do
* Do not delete any folder after downloading :)
* Do not put the target-region file in any format other than FASTA