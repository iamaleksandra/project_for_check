# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 13:41:13 2018

@author: aleksandra
"""

from setuptools import setup

setup(
    name='H_distribution_heliosphere',
    version='0.0.1',
    author='Volosatykh Aleksandra',
    description='Distribution function of neutral H within the heliosphere',
    py_modles=['mod1', 'statdistr2D'],
    # scripts=[''],
    test_suite='test',
    #install_requirements=['numpy>=1', 'matplotlib>=1']
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
    keywords='science astrophysics heliosphere'        
)

