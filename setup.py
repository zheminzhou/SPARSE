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
    version= '0.1.2',  
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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
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
    install_requires=['pycapnp', 'numpy', 'pandas', 'Cython', 'scipy', 'msgpack'],
    include_package_data=True,
    project_urls={ 
        'Bug Reports': 'https://github.com/zheminzhou/SPARSE/issues',
        'Source': 'https://github.com/zheminzhou/SPARSE/',
    },
)
