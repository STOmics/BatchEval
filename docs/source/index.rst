.. BatchEval documentation master file, created by
   sphinx-quickstart on Tue Sep  5 16:14:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BatchEval's documentation!
=====================================

BatchEval: Batch Effects Evaluation Workflow for Multi-batch Dataset Joint Analysis

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   tutorial/Tutorial


As genomic sequencing technology develops, multi-batch joint analysis of gene expression data can maximize the scientific value in the data set, supporting researchers in discovering more significant biological topics. However, joint analysis usually misses batch effects caused by abiotic deviations such as external noise in integrated data. When the batch effect makes an impact on the results of an experiment, the downstream analytical conclusion, if not the wrong experimental result, could have been incorrect. As a result, prior to data integration and processing, the batch effects in the data should be thoroughly investigated. We developed the BatchEval Pipeline, which is suitable for massive data integration, for evaluating batch effects in data and provide examination conclusions. The BatchEval Pipeline analyses multiple batches of data from multiple perspectives to detect how it affects of batch effects on integrated data. BatchEval Pipeline generates an HTML webpage report output for the assessment findings, which includes a basic explanation of the data, statistical analysis, a batch effect metric score, and visualization. BatchEval Pipeline analyses integrated data from multiple perspectives to determine whether it has to be corrected in an additional step employing batch effect removal approaches and how to do so, which is essential for downstream analysis.


