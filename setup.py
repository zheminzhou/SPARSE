"""A setuptools based setup module for SPARSE.

"""

from setuptools import setup, find_packages
from codecs import open
from os import path, walk, chdir, system

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='meta-sparse',  
    version= '0.1.6',  
    description='SPARSE indexes reference genomes in public databases into hierarchical clusters and uses it to predict origins of metagenomic reads.',
    long_description=long_description, 
    long_description_content_type='text/markdown',  
    url='https://github.com/zheminzhou/SPARSE/',
    author='Zhemin Zhou',  
    author_email= 'zhemin.zhou@warwick.ac.uk', 
    classifiers=[ 
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
#    py_modules=['SPARSE', 'parameter'],
    entry_points={
        'console_scripts': [
            'sparse = SPARSE.SPARSE:SPARSE',
        ],
    },
    keywords=['bioinformatics', 'microbial', 'metagenomics'], 
    package_data={'SPARSE':[
        'EM/*',
        'bin/*',
        'docs/*',
        'modules/*',
        'EM/Ipopt/include/*',
        'EM/Ipopt/lib64/*',
        'EM/Ipopt/*',
    ]},
    packages = ['SPARSE'], 
    package_dir = {'SPARSE':'.'},
    install_requires=['pycapnp==0.6.1', 'numpy==1.13.3', 'pandas==0.23.0', 'Cython', 'scipy==1.0.0', 'msgpack-python==0.4.8', 'msgpack==0.5.6'],
    include_package_data=True,
    project_urls={ 
        'Bug Reports': 'https://github.com/zheminzhou/SPARSE/issues',
        'Source': 'https://github.com/zheminzhou/SPARSE/',
    },
)
