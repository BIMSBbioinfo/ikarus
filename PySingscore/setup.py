from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='singscore',
      version='0.1',
      description='implement singscore (doi:10.1101/231217) in Python',
      url='http://github.com/kristyhoran/singscore',
      author='Kristy Horan',
      author_email='kristyhoran15@gmail.com',
      license='MIT',
      packages=['singscore', 'singscore.test'],
      install_requires=[
          'pandas','numpy','matplotlib','seaborn','scipy', 'statsmodels'
      ],
      zip_safe=False)