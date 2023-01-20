# scIBD
scIBD is a doublet detection tool for scCAS data.

scIBD is totally produced by Python.

The depending packages used in scIBD can be installed by the command pip/pip3 install -r requirements.txt 

Installation
-----

```bash
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
input: the count matrix(numpy or scipy)

output: the idx of predicted doublets; the doublet scores of all droplets

other parameters:

exprate: the expected calling rate of doublets, default 0.1




