.. BatchEval documentation master file, created by
   sphinx-quickstart on Tue Sep  5 16:14:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BatchEval's documentation!
=====================================

BatchEval: Batch Effects Evaluation Workflow for Multi-batch Dataset Joint Analysis

.. toctree::
   :maxdepth: 1

   installation
   tutorial/Tutorial

Overview
========
As genomic sequencing technology continues to advance, it becomes increasingly important to perform multiple dataset joint analysis of transcriptomics to understand complex biological systems. However, batch effect presents challenges for dataset integration, such as sequencing measured by different platforms and datasets collected at different times. Here, we develop a BatchEval Pipeline, which is used to evaluate batch effect of dataset integration and output a comprehensive report. This report consists of a series of HTML pages for the assessment findings, including a main page, a raw dataset evaluation page and several built-in methods evaluation pages. The main page exhibits basic information of integrated datasets, comprehensive score of batch effect and the most recommended method for batch effect removal to current datasets. The residual pages exhibit the evaluation details of raw dataset and evaluation results of many built-in batch effect removal methods after removing batch effect. This comprehensive report enables researchers to accurately identify and remove batch effect, resulting in more reliable and meaningful biological insights from integrated datasets. In summary, BatchEval Pipeline represents a significant advancement in batch effect evaluation and is a valuable tool to improve the accuracy and reliability of the experimental results.


