#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:35 PM
# @Author  : zhangchao
# @File    : trainer.py
# @Email   : zhangchao5@genomics.cn
import pandas as pd
import numpy as np
from typing import Union
from anndata import AnnData

from BatchEval.module.classifier import BatchClassifier


def domain_variance_score(merge_data: AnnData,
                          input_dims: int,
                          n_batch: int,
                          use_rep: str = "X_pca",
                          batch_key: str = "batch",
                          batch_size: int = 4096,
                          gpu: Union[str, int] = 0,
                          save_path: str = "./"):
    assert use_rep in merge_data.obsm_keys()
    classifier = BatchClassifier(input_dims=input_dims,
                                 n_batch=n_batch,
                                 data_x=merge_data.obsm[use_rep],
                                 batch_idx=merge_data.obs[batch_key].cat.codes,
                                 batch_size=batch_size,
                                 gpu=gpu)
    classifier.train(max_epochs=500, save_path=save_path)
    test_acc = classifier.test(pt_path=save_path)
    df = pd.DataFrame(data={"n_batch": n_batch,
                            "n_sample": merge_data.shape[0],
                            "Train Size": classifier.train_size,
                            "Accept Rate": np.around(1 - test_acc, decimals=4)},
                      index=["domain score"])
    return df
