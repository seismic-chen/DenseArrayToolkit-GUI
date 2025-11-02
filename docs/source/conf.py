import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

# 项目信息
project = 'DenseArrayToolkit-GUI'
copyright = '2024, Pengfei Zuo'
author = 'Pengfei Zuo, Yunfeng Chen, et al.'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
]

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']
