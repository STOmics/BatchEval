#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:23 PM
# @Author  : zhangchao
# @File    : kernel_plot.py
# @Email   : zhangchao5@genomics.cn
import matplotlib.pyplot as plt
import seaborn as sn
from io import BytesIO


def kernel_plot(merge_data, batch_key="batch", test_key="total_counts"):
    plt.figure(figsize=(12.8, 7.2))
    ax = plt.gca()
    ax.spines[["top", "right"]].set_visible(False)
    sn.kdeplot(x=test_key, data=merge_data.obs, hue=batch_key, common_norm=False)
    # plt.title(f"Kernel Fitting Curve by {test_key}".title())
    plt.xticks(fontsize=20, rotation=-35)
    plt.yticks(fontsize=20)
    plt.xlabel(test_key.title(), fontsize=20)
    plt.ylabel("Density", fontsize=20)
    plt.grid()
    fig_buffer = BytesIO()
    plt.savefig(fig_buffer, format="png", bbox_inches='tight', dpi=300)
    plt.close()
    return fig_buffer
