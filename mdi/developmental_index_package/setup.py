from setuptools import setup, find_packages

setup(
    name = 'developmental_index',
    version = '1.2',
    author = 'Ben Devlin',
    author_email = 'bad36@duke.edu',
    packages = find_packages('numpy', 'pandas', 'scipy'),
    url = 'https://github.com/bendevlin18/microglia-seq',
    scripts=["developmental_index.py"]

)