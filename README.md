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

if the input is an obeject of AnnData, the output is also an AnnData, the obs of returned AnnData adds two columns: obs['PredDBL'] is the predicted results where 1 indicates the predicted doublets and 0 indicates the singlets; obs['DBLscore'] is the doublet scores of all droplets.

if the input is the count matrix, the output are the idx of predicted doublets and the doublet scores of all droplets

-----
other parameters:

exprate: The expected calling rate of doublets, default 0.1.

strategy: The KNN graphing strategy, scIBD can adaptively opt a KNN graphing strategy. Users can also manually set it as 'PCoA' or 'PCA'.

core: The number of threads, default is the max core number depending on the terminals.

sim_rate: The ratio of simulated doublets in each iteration, default is 0.3.

nPC: The number of used principal components, default is 5.

neighbors: The number of neighbors used to construct KNN graph, default is 40.

n_tree: The number of trees in KNN constrcution, default is 30.





