#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:41 PM
# @Author  : zhangchao
# @File    : BatchEval.py
# @Email   : zhangchao5@genomics.cn
from collections import defaultdict

import scanpy as sc
import pandas as pd
import numpy as np
from typing import Union
from anndata import AnnData

from BatchEval.BatchEval import main_page, raw_page
from BatchEval.BatchEval.check_correct import check_correct
from BatchEval.BatchEval.check_raw import check_raw
from BatchEval.BatchEval.generate_pages import correct_page
from BatchEval.integration.bbknn import bbknn_integration
from BatchEval.integration.harmony import harmony_integration
from BatchEval.integration.spatiAlign import spatiAlign_integration
from BatchEval.utils import print_time, pca_lowrank

LINK = {"Harmony": "https://github.com/slowkow/harmonypy.git",
        "BBKNN": "https://github.com/Teichlab/bbknn.git",
        "spatiAlign": "https://github.com/STOmics/Spatialign.git"}


def _external(correct_data: AnnData,
              re_pca: bool = True,
              re_neigh: bool = True,
              use_rep: str = "X_pca",
              re_umap: bool = True):
    if re_pca:
        pca_lowrank(correct_data, use_rep=None, n_component=50)
    if re_neigh:
        sc.pp.neighbors(correct_data, use_rep=use_rep)
    if re_umap:
        sc.tl.umap(correct_data)


@print_time
def batch_eval(*data: AnnData,
               norm_log: bool = True,
               is_scale: bool = False,
               n_pcs: int = 50,
               n_neighbors: int = 50,
               batch_key: str = "batch",
               position_key: str = "X_umap",
               condition: Union[str, list, None] = None,
               count_key: str = "total_counts",
               celltype_key: Union[str, None] = None,
               external_list: list = None,
               report_path: str = "./",
               gpu: Union[str, int] = 0):
    """BatchEval Raw Dataset Pipeline

    Parameters
    -----------------
    *data: 'Anndata'
        Data matrix with rows for cells and columns for genes.
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
    condition: 'str, list, None'
        Label the experimental conditions. By default, the experimental conditions for each data are different.
    count_key: 'str'
    celltype_key: 'str'
    external_list: 'list'
        external integrated data list, and each element should be a dictionary, the dict key as follows,
        key:
            'correct_data',
            're_pca': whether to re-calculate PCA
            're_neigh': whether to re-calculate neighbors
            'use_rep': corrected embedding, and it should be saved in '*.obsm'
            're_umap': whether to re-calculate UMAP
            'batch_key'
    report_path: 'str'
    gpu: 'str', 'int'

    Return
    -----------------
    output_dict: 'dict'
    """
    total_score_df = pd.DataFrame(index=["BatchEval Score"])
    kbet_df = pd.DataFrame(index=["KBET Score", "95% p Value"])

    merge_data, data_dict_raw, img_dict_raw = check_raw(*data,
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

    raw_page(data_dict=data_dict_raw,
             img_dict=img_dict_raw,
             save_path=report_path,
             save_name="Raw Report.html")

    summary_score(data_dict_raw, name="Raw", df=total_score_df)
    summary_kbet(data_dict_raw, name="Raw", df=kbet_df)
    page_dict = {"Raw": "Raw report.html"}

    main_dict = defaultdict()
    if data_dict_raw["kbet"]["95% P Value"].values > 0.1:
        summary_describe = "The batch effect is much smaller than the biological variance " \
                           "and no additional batch effect processing is required."
    else:
        # spatiAlign
        spatialign_data = merge_data.copy()
        spatial_res = spatiAlign_integration(spatialign_data, batch_key, save_path=report_path, tau1=0.02, tau2=0.01,
                                             tau3=0.1)
        sc.pp.neighbors(spatial_res, use_rep="correct")
        sc.tl.umap(spatial_res)

        data_dict_spatialign, img_dict_spatialign = check_correct(spatial_res,
                                                                  batch_key="batch",
                                                                  use_rep="correct",
                                                                  position_key="X_umap",
                                                                  n_neighbors=n_neighbors,
                                                                  celltype_key=celltype_key,
                                                                  report_path=report_path,
                                                                  gpu=gpu)

        del spatial_res
        del spatialign_data
        summary_score(data_dict_spatialign, name="spatiAlign", df=total_score_df)
        summary_kbet(data_dict_spatialign, name="spatiAlign", df=kbet_df)
        correct_page(data_dict=data_dict_spatialign,
                     img_dict=img_dict_spatialign,
                     save_path=report_path,
                     save_name="spatiAlign Report.html",
                     correct_method="spatiAlign")

        page_dict.update({"spatiAlign": "spatiAlign Report.html"})

        # harmony
        harmony_data = merge_data.copy()
        harmony_integration(merge_data=harmony_data, batch_key="batch", n_pcs=n_pcs, use_rep="X_pca")
        sc.pp.neighbors(harmony_data, use_rep="X_pca_harmony")
        sc.tl.umap(harmony_data)
        data_dict_harmony, img_dict_harmony = check_correct(harmony_data,
                                                            batch_key="batch",
                                                            use_rep="X_pca_harmony",
                                                            position_key="X_umap",
                                                            n_neighbors=n_neighbors,
                                                            celltype_key=celltype_key,
                                                            report_path=report_path,
                                                            gpu=gpu)
        summary_score(data_dict_harmony, name="Harmony", df=total_score_df)
        summary_kbet(data_dict_harmony, name="Harmony", df=kbet_df)
        correct_page(data_dict=data_dict_harmony,
                     img_dict=img_dict_harmony,
                     save_path=report_path,
                     save_name="Harmony Report.html",
                     correct_method="Harmony")
        page_dict.update({"Harmony": "Harmony Report.html"})
        del harmony_data

        # bbknn
        bbknn_data = merge_data.copy()
        bbknn_integration(bbknn_data, batch_key="batch", n_pcs=n_pcs, use_rep="X_pca")
        sc.tl.umap(bbknn_data)
        data_dict_bbknn, img_dict_bbknn = check_correct(bbknn_data,
                                                        batch_key="batch",
                                                        use_rep="X_pca",
                                                        position_key="X_umap",
                                                        n_neighbors=n_neighbors,
                                                        celltype_key=celltype_key,
                                                        report_path=report_path,
                                                        gpu=gpu)

        correct_page(data_dict=data_dict_bbknn,
                     img_dict=img_dict_bbknn,
                     save_path=report_path,
                     save_name="BBKNN Report.html",
                     correct_method="BBKNN")

        page_dict.update({"BBKNN": "BBKNN Report.html"})

        summary_score(data_dict_bbknn, name="BBKNN", df=total_score_df)
        summary_kbet(data_dict_bbknn, name="BBKNN", df=kbet_df)
        del bbknn_data

        if external_list is not None:
            for idx, ext in enumerate(external_list):
                assert isinstance(ext, dict)
                check_k = []
                for k in ["correct_data", "re_pca", "re_neigh", "use_rep", "re_umap", "batch_key"]:
                    if k not in ext.keys():
                        check_k.append(k)
                if len(check_k) != 0:
                    raise ValueError(f"`external_list` got invalid values: {check_k}, please checkout.")
                _external(correct_data=ext["correct_data"],
                          re_pca=ext["re_pca"],
                          re_neigh=ext["re_neigh"],
                          re_umap=ext["re_umap"],
                          use_rep=ext["use_rep"])

                data_dict_ext, img_dict_ext = check_correct(ext["correct_data"],
                                                             batch_key=ext["batch_key"],
                                                             use_rep=ext["use_rep"],
                                                             position_key="X_umap",
                                                             n_neighbors=n_neighbors,
                                                             celltype_key=celltype_key,
                                                             report_path=report_path,
                                                             gpu=gpu)
                correct_page(data_dict=data_dict_ext,
                             img_dict=img_dict_ext,
                             save_path=report_path,
                             save_name=f"External{idx} Report.html",
                             correct_method=f"External{idx}")

                summary_score(data_dict_bbknn, name=f"External{idx}", df=total_score_df)
                summary_kbet(data_dict_bbknn, name=f"External{idx}", df=kbet_df)
                page_dict.update({f"External{idx}": f"External{idx} Report.html"})

        tmp_method = total_score_df.columns[total_score_df.loc['BatchEval Score'].argmax()]
        summary_describe = f"This data has batch effect and requires further processing, " \
                           f"and recommend '{tmp_method}'. More details of '{tmp_method}' can be found in '{LINK[tmp_method]}'."
        main_dict.update({"summary-kbet": kbet_df,
                          "summary-score": total_score_df})

    main_dict.update({"summary-conclusion": pd.DataFrame(columns=[f"{summary_describe}"])})
    main_dict.update({"describe": data_dict_raw["describe"]})
    main_page(
        data_dict=main_dict,
        img_dict=img_dict_raw,
        save_path=report_path,
        save_name="BatchEval Report.html",
        pages_dict=page_dict)


def summary_score(data_dict, name, df):
    df[name] = (data_dict["domain"]["Accept Rate"].values
                + data_dict["lisi"]["LISI(f1 score)"].values
                + data_dict["silhouette"]["SS_f1"].values) / 3


def summary_kbet(data_dict, name, df):
    df[name] = np.array([data_dict["kbet"]["Accept Rate"].values, data_dict["kbet"]["95% P Value"].values]).flatten()
