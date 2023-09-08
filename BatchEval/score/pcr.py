#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/8 17:09
# @Author  : zhangchao
# @File    : pcr.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import LinearRegression


def pcr_comparison(adata_pre, adata_post, covariate, embed=None, n_comps=50, scale=True, verbose=False):
    """Compare the explained variance before and after integration.

    Return either the difference of variance contribution before and after integration or a score between 0 and 1
    (``scaled=True``) with 0 if the variance contribution hasn't changed. The larger the score,
    the more different the variance contributions are before and after integration.

    Parameters
    ----------
    adata_pre: anndata object before integration
    adata_post: anndata object after integration
    covariate: Key for ``adata_post.obs`` column to regress against
    embed
    n_comps
    scale
    verbose

    Returns
    -------

    """
    if embed == "X_pca":
        embed = None

    pcr_before = pcr(
        adata_pre,
        covariate=covariate,
        recompute_pca=True,
        n_comps=n_comps,
        verbose=verbose,
    )

    pcr_after = pcr(
        adata_post,
        covariate=covariate,
        embed=embed,
        recompute_pca=True,
        n_comps=n_comps,
        verbose=verbose,
    )

    if scale:
        score = (pcr_before - pcr_after) / pcr_before
        if score < 0:
            print(
                "Variance contribution increased after integration!\n"
                "Setting PCR comparison score to 0."
            )
            score = 0
        return score
    else:
        return pcr_after - pcr_before


def pcr(adata, covariate, embed=None, n_comps=50, recompute_pca=True, verbose=False):

    if verbose:
        print(f"covariate: {covariate}")
    covariate_values = adata.obs[covariate]

    # use embedding for PCA
    if embed is not None:
        assert embed in adata.obsm
        if verbose:
            print(f"Compute PCR on embedding n_comps: {n_comps}")
        return pc_regression(adata.obsm[embed], covariate_values, n_comps=n_comps)

    # use existing PCA computation
    elif (recompute_pca is False) and ("X_pca" in adata.obsm) and ("pca" in adata.uns):
        if verbose:
            print("using existing PCA")
        return pc_regression(
            adata.obsm["X_pca"], covariate_values, pca_var=adata.uns["pca"]["variance"]
        )

    # recompute PCA
    else:
        if verbose:
            print(f"compute PCA n_comps: {n_comps}")
        return pc_regression(adata.X, covariate_values, n_comps=n_comps)


def pc_regression(data, covariate, pca_var=None, n_comps=50, svd_solver="arpack", verbose=False):
    if isinstance(data, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix)):
        matrix = data
    else:
        raise TypeError(
            f"invalid type: {data.__class__} is not a numpy array or sparse matrix"
        )

    # perform PCA if no variance contributions are given
    if pca_var is None:

        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = "full"
            # convert to dense bc 'full' is not available for sparse matrices
            if sparse.issparse(matrix):
                matrix = matrix.todense()

        if verbose:
            print("compute PCA")
        X_pca, _, _, pca_var = sc.tl.pca(
            matrix,
            n_comps=n_comps,
            use_highly_variable=False,
            return_info=True,
            svd_solver=svd_solver,
            copy=True,
        )
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]

    # PC Regression
    if verbose:
        print("fit regression on PCs")

    # handle categorical values
    if pd.api.types.is_numeric_dtype(covariate):
        covariate = np.array(covariate).reshape(-1, 1)
    else:
        if verbose:
            print("one-hot encode categorical values")
        covariate = pd.get_dummies(covariate)

    # fit linear model for n_comps PCs
    r2 = []
    for i in range(n_comps):
        pc = X_pca[:, [i]]
        lm = LinearRegression()
        lm.fit(covariate, pc)
        r2_score = np.maximum(0, lm.score(covariate, pc))
        r2.append(r2_score)

    Var = pca_var / sum(pca_var) * 100
    R2Var = sum(r2 * Var) / 100

    return R2Var
