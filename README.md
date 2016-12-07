# potato_pigs_analyses

# Installation

To get a conda environment setup run
```
conda create -n ppps pip python=3 numpy
source activate ppps
conda install pandas scipy matplotlib seaborn h5py
pip install biom-format
```

To install the packages required for emperor run

`conda install -c conda-forge emperor h5py`

To install the packages required for the balance trees run 

`pip install git+https://github.com/biocore/gneiss.git`

To install canvas used for the permutative ANCOM tests run

`pip install git+https://github.com/mortonjt/canvas.git`
