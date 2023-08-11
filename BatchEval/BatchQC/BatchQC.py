#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:41 PM
# @Author  : zhangchao
# @File    : BatchQC.py
# @Email   : zhangchao5@genomics.cn
from typing import Union
from anndata import AnnData

from BatchEval.BatchQC.batchqc_raw import batchqc_raw
from BatchEval.utils import print_time


@print_time
def batch_qc(*data: AnnData,
             qc_mode: str = "raw",
             adjust_method: str = "harmony",
             norm_log: bool = True,
             is_scale: bool = False,
             n_pcs: int = 50,
             n_neighbors: int = 50,
             batch_key: str = "batch",
             position_key: str = "X_umap",
             use_rep: str = "X_pca",
             condition: Union[str, list, None] = None,
             count_key: str = "total_counts",
             celltype_key: Union[str, None] = None,
             report_path: str = "./",
             gpu: Union[str, int] = 0):
    """BatchQC Raw Dataset Pipeline

    Parameters
    -----------------
    *data: 'Anndata'
        Data matrix with rows for cells and columns for genes.
    qc_mode: 'str'
        select QC mode. optional ['raw', 'adjust']
    adjust_method: 'str'
        select batch effects removal approach. It works for 'qc_mode = adjust'.
    norm_log: 'bool'
        Whether to preprocess data. 'sc.pp.normalization()', 'sc.pp.log1p()'
    is_scale: 'bool'
         Whether to preprocess data. 'sc.pp.scale()'
    n_pcs: 'int'
        Number of principal components retained in PCA. default, 50.
    n_neighbors: 'int'
        Calculate the nearest neighbors of a local area. default, 100.
    batch_key: 'str'
        Label the data batches.
    position_key: 'str'
        Compute the coordinate space of the nearest neighbor.
    use_rep: 'str'
        the embedding be used to batch effect removal, default, 'X_pca'. It works for 'qc_mode = adjust'.
    condition: 'str, list, None'
        Label the experimental conditions. By default, the experimental conditions for each data are different.
    count_key: 'str'
    celltype_key: 'str'
    report_path: 'str'
    gpu: 'str', 'int'

    Return
    -----------------
    output_dict: 'dict'
    """
    if qc_mode == "raw":
        batchqc_raw(
            *data,
            norm_log=norm_log,
            is_scale=is_scale,
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            batch_key=batch_key,
            position_key=position_key,
            condition=condition,
            count_key=count_key,
            celltype_key=celltype_key,
            report_path=report_path,
            gpu=gpu)
    elif qc_mode == "adjust":
        raise ValueError("This particular feature is currently undergoing development and refinement...")
    else:
        raise ValueError("Got an invalid mode, which only support 'raw' and 'adjust'!")
