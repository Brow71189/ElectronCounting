#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 14:55:28 2017

@author: mittelberger2
"""

from distutils.core import setup
from distutils.extension import Extension

setup(
      name = 'electron counting',
      py_modules = ['ElectronCounting.electron_counting'],
      ext_modules = [Extension('ElectronCounting.c_electron_counting', ['ElectronCounting/c_electron_counting.c'])],
)
