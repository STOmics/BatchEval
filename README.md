[![python >3.8.8](https://img.shields.io/badge/python-3.8.8-brightgreen)](https://www.python.org/)

# BatchEval: Batch Effects Evaluation Workflow for Multi-batch Dataset Joint Analysis

As genomic sequencing technology develops, multi-batch joint analysis of gene expression data can maximize the
scientific value in the data set, supporting researchers in discovering more significant biological topics. However,
joint analysis usually misses batch effects caused by abiotic deviations such as external noise in integrated data. When
the batch effect makes an impact on the results of an experiment, the downstream analytical conclusion, if not the wrong
experimental result, could have been incorrect. As a result, prior to data integration and processing, the batch effects
in the data should be thoroughly investigated. We developed the BatchEval Pipeline, which is suitable for massive data
integration, for evaluating batch effects in data and provide examination conclusions. The BatchEval Pipeline analyses
multiple batches of data from multiple perspectives to detect how it affects of batch effects on integrated data.
BatchEval Pipeline generates an HTML webpage report output for the assessment findings, which includes a basic
explanation of the data, statistical analysis, a batch effect metric score, and visualization. BatchEval Pipeline
analyses integrated data from multiple perspectives to determine whether it has to be corrected in an additional step
employing batch effect removal approaches and how to do so, which is essential for downstream analysis.

# Dependences

[![torch-1.10.0](https://img.shields.io/badge/torch-1.10.0-red)](https://pytorch.org/get-started/previous-versions/)
[![pandas-1.2.4](https://img.shields.io/badge/pandas-1.2.4-lightgrey)](https://github.com/pandas-dev/pandas)
[![scikit-learn-0.24](https://img.shields.io/badge/scikit-0.24.x-brightgreen)](https://github.com/scikit-learn/scikit-learn/tree/0.24.X)
[![scipy-0.12.x](https://img.shields.io/badge/scipy-0.12.x-yellow)](https://github.com/scipy/scipy/tree/maintenance/0.12.x)
[![distinctipy-1.2.2](https://img.shields.io/badge/distinctipy-1.2.2-green)](https://github.com/alan-turing-institute/distinctipy/tree/v1.2.2)
[![lxml-4.9.2](https://img.shields.io/badge/lxml-4.9.2-9cf)](https://github.com/lxml/lxml/tree/lxml-4.9.2)                  
[![scanpy-1.9.1](https://img.shields.io/badge/scanpy-1.9.1-informational)](https://pypi.org/project/scanpy/)

# Install
   
```python
pip install BatchEval
```    
        
```git
git clone https://github.com/STOmics/BatchEval.git

cd BatchEval

python setup.py install
```

# Tutorial

- BatchEval Pipeline           
  Check the batch effects of raw data.
  ```python
  import os
  import scanpy as sc
  
  from BatchEval import batch_eval
  ```
  
  ```python
  data = [sc.read_h5ad(os.path.join("./demo_data", t)) for t in sorted(os.listdir("./demo_data")) if t.endswith("h5ad")]
  ```
        
  ```python
  batch_eval(*data,
          qc_mode="raw",
          norm_log=True,
          is_scale=True,
          n_pcs=50,
          n_neighbors=100,
          batch_key="batch",
          position_key="X_umap",
          condition=None,
          use_rep="X_pca",
          count_key="total_counts",
          celltype_key=None,
          report_path="./output")
  ```      
        
- **Note: 1) the `data*` is `AnnData` format. 2) In the batchQC pipeline, the original expression will be counted, so
  the input expression matrix should be the original captured expression.**        
  Check the batch effects of adjust data.

# Disclaimer

***This is not an official product.***       
         
        


            
            
            
            
