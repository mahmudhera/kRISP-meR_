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
```buildoutcfg
sudo apt-get update -y
sudo apt-get install -y samtools
```
* **Biopython**: install using: `pip install biopython`
* **Python binding of Jellyfish** (gmarics project). This is quite tricky. Need to assess this in detail later. Tried to do the following:
```buildoutcfg
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
```
Python binding means that if you try to import jellyfish from a python script (the code is: `import jellyfish`), that will work and you will be able to invoke Jellyfish program from python.
If the installation steps mentioned above does not bind with python, (although the documentation of Jellyfish does say that this should): you may try the swig binding instructions from https://github.com/gmarcais/Jellyfish/blob/master/swig/Readme.md).
```buildoutcfg
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
1. **Download** the github directory
1. **Locate** your sequenced reads file (can be FASTA or FASTQ) and the target-region file (_must_ be FASTA); these files can be anywhere on the filesystem
1. **Run** the tool with the command: `python krispmer.py <reads_filename> <target_filename> <num_mismatches> <output_file>`
The arguments, alongwith other options can be seen using `python krispmer.py -h`.

## Do not do
1. Do not delete any folder after downloading :)
1. Do not put the target-region file in any format other than FASTA

## Options available in kRISP-mER
