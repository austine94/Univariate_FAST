Readme for Sequential Detection of Emergent Anomalous Structures in Functional Data code

The zip file contains three things.
	1) All code needed to recreate each of the results in the paper (using main.R).
	2) A specific function, Build_Table1.R, to satisfy the reproducibility requirement.
	3) A specific function, Build_Table6.R, to reproduce one of the comparative studies (as requested in an editorial comment).
Running each of these three elements will be described in turn below.

0) Package Requirements:
	These functions will require that the R packages "fda", "fdaoutlier", "aplpack", "tidyverse", and "mvtnorm" are installed.
	Additionally, storing the tables will require the R package "gridExtra", however the results could be printed to console without.

1) All Results:
	Set the working directory to the root of the unzipped folder.
	Open main.R
	Running main.R will produce .csv files of every table in the paper in the Results folder.
2) Build Table 1
	Set the working directory to the root of the unzipped folder.
	Open Build_Table1.R
	Running this script will save a .pdf version of Table 1 in the Results folder.
3) Build Table 6
	Set the working directory to the root of the unzipped folder.
	Open Build_Table6.R
	Running this script will save a .pdf version of Table 13 in the Results folder.

