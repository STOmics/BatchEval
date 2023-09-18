[![python >3.8.8](https://img.shields.io/badge/python-3.8.8-brightgreen)](https://www.python.org/)
[![Downloads](https://static.pepy.tech/badge/BatchEval)](https://pepy.tech/project/BatchEval)
[![Documentation Status](https://readthedocs.org/projects/batcheval/badge/?version=latest)](https://batcheval.readthedocs.io/en/latest/?badge=latest)
# BatchEval: Batch Effects Evaluation Workflow for Multi-batch Dataset Joint Analysis

As genomic sequencing technology continues to advance, it becomes increasingly important to perform multiple dataset joint analysis of transcriptomics to understand of complex biological systems. However, batch effect removal present challenges for data integration, such as tissue sequencing measured by different platforms, collected at different times. Here, we developed a BatchEval Pipeline, which is used to evaluate batch effect of data integration and output a comprehensive summarized of batch effect. It consists of a series HTML page for the assessment findings, including a main page, raw dataset evaluation page and a series built-in methods of data integration evaluation page. The main page exhibition basic information of integrated data, comprehensive score of batch effect and recommend the best method for batch effect removal of current data.  The residual pages are exhibition the evaluation details of raw dataset and using built-in methods after remove batch effect, respectively. This comprehensive report enables researchers to accurately identify and remove batch effect, resulting in more reliable and meaningful biological insights from integrated datasets. In summary, BatchEval Pipeline represents a significant advancement in the field of batch effect evaluation and provides a valuable tool for researchers to improve the accuracy and reliability of their experimental results.

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
or        
```git
git clone https://github.com/STOmics/BatchEval.git

cd BatchEval

python setup.py install
```

# Tutorial
Quick Start [https://batcheval.readthedocs.io/en/latest/?badge=latest](https://batcheval.readthedocs.io/en/latest/?badge=latest)

# Disclaimer

***This is not an official product.***       
         
        


            
            
            
            
