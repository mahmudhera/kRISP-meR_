# kRISP-mER
Reference free guide RNA designing tool

# Restrictions
1. Must keep the personalized, inputs and output folders
2. Must give the reads and targets in the inputs directory

# Requirements
1. sklearn: pip install scikit-learn==0.16.1
2. Biopython: pip install biopython
3. Python binding of Jellyfish (gmarics project). This is quite tricky. Need to assess this in detail later. Tried to do the following:
```buildoutcfg
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
```
Python binding means that if you try to import jellyfish from a python script (the code is: `import jellyfish`), that will work and you will be able to invoke Jellyfish program from python.

The installation steps mentioned above does not bind with python for some reason (although the documentation of Jellyfish does say that this should. Therefore, had to try the swig binding instructions from https://github.com/gmarcais/Jellyfish/blob/master/swig/Readme.md) This still does not seem to solve the issue. The code that I ran is as follows:

```buildoutcfg
cd swig/python
python setup.py build
python setup.py install --prefix=$PREFIX
```

Python binding issue still remains. Will resolve this later.
