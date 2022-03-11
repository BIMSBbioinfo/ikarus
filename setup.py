from setuptools import setup

with open("readme.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='ikarus',
    version='0.0.2',
    packages=['ikarus'],
    desription='Machine Learning classifier of tumor cells',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'scanpy',
        'anndata',
        'pyscenic'
    ]
)
