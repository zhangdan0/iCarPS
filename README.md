# iCarPS
iCarPS provides a machine-learning method for carbonylation site prediction. It is implemented by WEKA library--random forest. At present, it provides user-friendly webserver to predict carbonylation sites freely, refering to http://lin-group.cn/server/iCarPS/webServer.html. Here we also provide our source codes for researchers who would like to run locally and conveniently to predict the query protein sequences in bulk.

# Installation
Download iCarPS by 

`https://github.com/zhangdan0/iCarPS.git`

# Requirements
* Python>=3.0
* Java

Since the package is written in python 3.7, please use version 3.0 or above of python and the pip tool must be installed first. iCarPS uses the following dependencies: os, sys, re, numpy, pandas. You can install these packages first, by the following commands:
```
pip install os 
pip install sys
pip install re
pip install pandas 
pip install numpy
```
Besides, you must ensure that your local computer has a Java environment before running. You can open the cmd command window and input to check it:

`java -version`

# Useage
* For users who want to perform carbonylation site prediction by our provided model :
* cd to the ./iCarPS_offline folder which contains iCarPS_offline.py and run as:

`python iCarPS_offline.py [PredictType] [query sequences file]`

1. the first parameter [PredictType] is used to assign a suitable type for prediction, including K, P, R, or T.
2. the second parameter [query sequences file] is the input sequence file, which must be fasta format. Please refer to those files in ./input/example folder.

## Example:

`python iCarPS_offline.py K ./input/test.txt`

Finally, the result file is in output folder, named '*_finalresult.txt'
