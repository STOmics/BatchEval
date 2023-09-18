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
   api/

Overview
========
As genomic sequencing technology continues to advance, it becomes increasingly important to perform multiple dataset joint analysis of transcriptomics to understand of complex biological systems. However, batch effect removal present challenges for data integration, such as tissue sequencing measured by different platforms, collected at different times. Here, we developed a BatchEval Pipeline, which is used to evaluate batch effect of data integration and output a comprehensive summarized of batch effect. It consists of a series HTML page for the assessment findings, including a main page, raw dataset evaluation page and a series built-in methods of data integration evaluation page. The main page exhibition basic information of integrated data, comprehensive score of batch effect and recommend the best method for batch effect removal of current data.  The residual pages are exhibition the evaluation details of raw dataset and using built-in methods after remove batch effect, respectively. This comprehensive report enables researchers to accurately identify and remove batch effect, resulting in more reliable and meaningful biological insights from integrated datasets. In summary, BatchEval Pipeline represents a significant advancement in the field of batch effect evaluation and provides a valuable tool for researchers to improve the accuracy and reliability of their experimental results.


