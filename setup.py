from setuptools import setup

setup(
    name='ikarus',
    version='0.0.1',
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
