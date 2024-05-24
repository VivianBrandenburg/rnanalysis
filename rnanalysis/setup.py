#!python
from setuptools import setup, find_packages

setup(
    name='rnaseq_analysis',
    packages=find_packages(),
    version='0.1.0',
    author="VivianBrandenburg",
    install_requires=[
        'pandas', 'numpy', 'matplotlib', 'seaborn']
)
