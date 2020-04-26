# Expected-Genetic-Similarity-Matrices
Summary
-------

<!-- [![Build Status](https://travis-ci.org/alexandrebouchard/phylosmcsampler.png?branch=master)](https://travis-ci.org/alexandrebouchard/phylosmcsampler) -->

Expected genetic similarity matrices is a recently proposed approach for computing genetic similarities using phylogenies. 
The standard method for computing similarity matrices involves the  
inner product of observed genetic variant matrices. Such an approach is inaccurate if genotypes are not available, or not densely sampled, or of poor quality (for example, genetic analysis of extinct species). We provide a new method for computing genetic similarities among individuals using phylogenetic trees. Our method can supplement (or stand in for) computations based on genetic sequences. 

Installation
------------


There are several options available to install the package:

### Compile using the provided gradle script

- Check out the source ``git clone https://github.com/shijiaw/Expected-Genetic-Similarity-Matrices.git``
- Compile using ``./gradlew installDist``
- Add the jars in ``build/install/Expected-Genetic-Similarity-Matrices/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone https://github.com/shijiaw/Expected-Genetic-Similarity-Matrices.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Usage
-----

### Quick start
Set the directory of your input tree in src/GSM/gsm.java

-String inputTree = " ";

Set the directory of your output file for EGSM in src/GSM/gsm.java

-String OutputFileName = " ";

Note that the input tree is in .newick format and the names of individuals are of form ``leaf_1, leaf_2,...''.

Run src/GSM/gsm.java to compute the EGSM. 

