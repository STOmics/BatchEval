#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:11 PM
# @Author  : zhangchao
# @File    : setup.py.py
# @Email   : zhangchao5@genomics.cn
import setuptools
from wheel.bdist_wheel import bdist_wheel

import BatchEval

__version__ = BatchEval.__version__


class BDistWheel(bdist_wheel):
    def get_tag(self):
        return (self.python_tag, "none", "any")


cmdclass = {
    "bdist_wheel": BDistWheel,
}

requirements = open("requirements.txt").readline()

setuptools.setup(
    name="BatchEval",
    version=__version__,
    author="zhangchao",
    author_email="zhangchao5@genomics.cn",
    url="https://github.com/STOmics/BatchEval.git",
    description="BatchEval: Batch Effects Evaluation Workflow for Multi-batch Dataset Joint Analysis",
    python_requires=">=3.8",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    cmdclass=cmdclass,
    package_data={'': ["*.html"]},
)
