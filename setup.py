from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import os

this_directory = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(this_directory, 'README.md')):
    with open(os.path.join(this_directory, 'README.md'), 'r') as f:
        long_description = f.read()
else:
    long_description = '''## Arthorian quest\n'''

if os.path.exists(os.path.join(this_directory, 'requirements.txt')):
    with open(os.path.join(this_directory, 'requirements.txt'), 'r') as f:
        requirements = [line.split('#')[0].strip() for line in f.readlines()]
        requirements = [line for line in requirements if line]
else:
    requirements = []

setup(
    name='arthorian-quest',
    version='0.1.7',
    description='using Arthor and filtering the results with Fragmenstein',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    #extras_require={'jupyter': ['jupyter']},
    url='https://github.com/matteoferla/arthorian-quest',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    classifiers=[ # https://pypi.org/classifiers/
        'Development Status :: 3 - Alpha', # Development Status :: 5 - Production/Stable
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)