#!/usr/bin/env python3

from setuptools import setup

config = {
    'name': 'ppcg-qc-from-sanger',
    'description': 'A tool to extract PPCG defined QC metrics from Sanger pipeline results',
    'author': 'Yaobo Xu',
    'url': 'https://github.com/cancerit/ppcg-qc-from-sanger',
    'download_url': '',
    'author_email': 'cgphelp@sanger.ac.uk',
    'version': '0.4.1',
    'python_requires': '>= 3.6',
    'setup_requires': ['pytest'],
    'install_requires': [],
    'packages': ['ppcg_qc_from_sanger'],
    'package_data': {'ppcg_qc_from_sanger': ['config/*.json']},
    'entry_points': {
        'console_scripts': ['ppcg-qc-from-sanger=ppcg_qc_from_sanger.command_line:main'],
    },
}

setup(**config)
