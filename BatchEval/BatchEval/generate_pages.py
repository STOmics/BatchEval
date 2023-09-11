#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/7 14:04
# @Author  : zhangchao
# @File    : generate_pages.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import os
import os.path as osp
import pwd
import time
import pkgutil
import pandas as pd
from lxml import etree

from BatchEval.utils import embed_text, embed_tabel
from BatchEval.utils.html_utils import embed_link_table, embed_table_imgs, embed_single_link


def main_page(data_dict: dict,
              save_path: str,
              save_name: str = "BatchEval Report.html",
              pages_dict: dict = {"Raw": "Raw_report.html"}):
    """Generate HTML report Main page

    Parameters
    ----------
    data_dict: display dataset
    save_path:
    save_name:
    pages_dict: subpage dict, include 'raw page', 'evaluation page' ...

    Returns
    -------
    """

    html = etree.HTML(pkgutil.get_data("BatchEval", "report_template_main.html").decode())
    embed_basic(html)

    allow_keys = ["describe", "summary-kbet", "summary-score", "summary-conclusion"]
    for key in allow_keys:
        assert key in data_dict.keys()

    assert isinstance(data_dict["describe"], pd.DataFrame)
    embed_tabel(data_dict["describe"], html, pos="h4", name="describe", is_round=False)
    embed_tabel(data_dict["summary-kbet"], html, pos="div", name="div-kbet", is_round=True)
    embed_tabel(data_dict["summary-score"], html, pos="div", name="div-score", is_round=True)
    embed_tabel(data_dict["summary-conclusion"], html, pos="div", name="div-conclusion", is_round=False)

    embed_link_table(html, pos="div", name="div-content", link_dict=pages_dict)

    tree = etree.ElementTree(html)
    os.makedirs(save_path, exist_ok=True)
    tree.write(osp.join(save_path, f"{save_name}"))
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} The '{save_name}' has been saved to {save_path}")


def raw_page(data_dict: dict,
             img_dict: dict,
             save_path: str,
             save_name: str = "Raw_report.html"):
    """Generate HTML report subpage

    Parameters
    ----------
    data_dict
    img_dict
    save_path
    save_name

    Returns
    -------

    """
    html = etree.HTML(pkgutil.get_data("BatchEval", "report_template_raw.html").decode())
    embed_basic(html)
    embed_statistical_table(data_dict, html)
    embed_biological_table(data_dict, html)

    embed_statistical_img(img_dict, html)
    embed_biological_img(img_dict, html)

    embed_single_link(tree=html, pos="button", name="back2main", link_dict={"Go Back...": "BatchEval Report.html"})

    tree = etree.ElementTree(html)
    os.makedirs(save_path, exist_ok=True)
    tree.write(osp.join(save_path, f"{save_name}"))
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} The '{save_name}' has been saved to {save_path}")


def correct_page(
        data_dict: dict,
        img_dict: dict,
        save_path: str,
        save_name: str,
        correct_method: str = None):
    html = etree.HTML(pkgutil.get_data("BatchEval", "report_template_adjust.html").decode())
    embed_basic(html)

    if correct_method is not None:
        embed_text(html, pos="h1", name="name", text=f"BatchEval Report ({correct_method})")

    embed_biological_table(data_dict, html)

    embed_biological_img(img_dict, html)

    embed_single_link(tree=html, pos="button", name="back2main", link_dict={"Go Back...": "BatchEval Report.html"})

    tree = etree.ElementTree(html)
    os.makedirs(save_path, exist_ok=True)
    tree.write(osp.join(save_path, f"{save_name}"))
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} The '{save_name}' has been saved to {save_path}")


def embed_statistical_img(img_dict, html):
    allowed_keys2 = ["box", "cdf", "kernel", "var-mean"]
    tmp_key = []
    for key in allowed_keys2:
        if key in img_dict.keys():
            continue
        else:
            tmp_key.append(key)
    if len(tmp_key) != 0:
        print(f"{tmp_key} can not found in `data_dict.keys()` [{img_dict.keys()}]")
        raise ValueError

    src_dict = {
        "Kernel Distribution Curve of UMICount Total": img_dict["kernel"],
        "CDF Curve of UMICount Total": img_dict["cdf"]
    }
    embed_table_imgs(src_dict, tree=html, pos="div", class_name="cal-curve")

    src_dict = {
        "Box-plot of UMICount Total": img_dict["box"],
        "UMICount Scatter Plot for Each Spot": img_dict["var-mean"]
    }
    embed_table_imgs(src_dict, tree=html, pos="div", class_name="cal-plot")


def embed_biological_img(img_dict, html):
    allowed_keys2 = ["umap_batch", "umap_type", "joint"]
    tmp_key = []
    for key in allowed_keys2:
        if key in img_dict.keys():
            continue
        else:
            tmp_key.append(key)
    if len(tmp_key) != 0:
        print(f"{tmp_key} can not found in `data_dict.keys()` [{img_dict.keys()}]")
        raise ValueError

    src_dict = {
        "Joint": img_dict["joint"],
        "UMAP-Batch": img_dict["umap_batch"],
        "UMAP-Type": img_dict["umap_type"]
    }
    embed_table_imgs(buffer_dict=src_dict, tree=html, pos="div", class_name="bio-plot")


def embed_statistical_table(data_dict, html):
    allowed_keys1 = ["f-test", "ks-test", "confound"]

    tmp_key = []
    for key in allowed_keys1:
        if key in data_dict.keys():
            continue
        else:
            tmp_key.append(key)
    if len(tmp_key) != 0:
        print(f"{tmp_key} can not found in `data_dict.keys()` [{data_dict.keys()}]")
        raise ValueError

    embed_tabel(data_dict["f-test"], html, pos="h4", name="f-test")
    embed_tabel(data_dict["ks-test"], html, pos="h4", name="ks-test")
    embed_tabel(data_dict["confound"], html, pos="h4", name="confound")


def embed_biological_table(data_dict, html):
    allowed_keys1 = ["domain", "kbet", "lisi", "silhouette"]

    tmp_key = []
    for key in allowed_keys1:
        if key in data_dict.keys():
            continue
        else:
            tmp_key.append(key)
    if len(tmp_key) != 0:
        print(f"{tmp_key} can not found in `data_dict.keys()` [{data_dict.keys()}]")
        raise ValueError

    embed_tabel(data_dict["domain"], html, pos="h4", name="domain")
    embed_tabel(data_dict['kbet'], html, pos="h4", name="kbet")
    embed_tabel(data_dict['lisi'], html, pos="h4", name="lisi")
    embed_tabel(data_dict['silhouette'], html, pos="h4", name="silhouette")


def embed_basic(html):
    embed_text(html,
               pos="h4",
               name="username",
               text=f"Report By: {pwd.getpwuid(os.getuid())[0]}")

    embed_text(html,
               pos="h5",
               name="runtime",
               text=f"Report Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
