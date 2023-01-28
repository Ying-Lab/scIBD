# scIBD
scIBD is a doublet detection tool for scCAS data.

scIBD is totally produced by Python.

The depending packages used in scIBD can be installed by the command pip/pip3 install -r requirements.txt 

Installation
-----

```bash
pip install -r requirements.txt 
git clone git://github.com/Ying-Lab/scIBD
cd scIBD
python setup.py install

```
Running
-----
```bash
import scibd as si
KNNITER = si.KNNIter(input)
result = KNNITER.IterCall()

```
Parameters
-----
input: the AnnData; or the count matrix(numpy or scipy)

output: 

if the input is an obeject of AnnData, the output is also an AnnData, the obs adds two columns: obs['PredDBL'] is the predicted results where 1 indicates the predicted doublets and 0 indicates the singlets; obs['DBLscore'] is the doublet scores of all droplets.

if the input is the count matrix, the output are the idx of predicted doublets and the doublet scores of all droplets

-----
other parameters:

exprate: the expected calling rate of doublets, default 0.1




