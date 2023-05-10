from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='ikarus',
    version='0.0.3',
    packages=['ikarus'],
    desription='Machine Learning classifier of tumor cells',
    license='MIT',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'scanpy',
        'anndata',
        'pyscenic',
        'ctxcore'
    ]
)
