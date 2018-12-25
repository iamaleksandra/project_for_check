# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name='H_distribution_heliosphere',
    version='0.0.1',
    author='Volosatykh Aleksandra',
    description='Distribution function of neutral H within the heliosphere',
    py_modles=['mod1', 'statdistr2D'],
    test_suite='test',
    install_requirements=['numpy>=1.9', 'matplotlib>=3.0.2'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords='science astrophysics heliosphere'        
)

