#!/usr/bin/env python

from distutils.core import setup

import pyrexplorer

setup(
    name='pyReXplorer',
    version=pyrexplorer.__version__,
    author=pyrexplorer.__author__,
    author_email=pyrexplorer.__contact__,
    packages=['pyrexplorer'],
    license='Apache License, Version 2.0',
    description='Python package aimed to explore relations in analyzed data',
    long_description=open('README').read(),
)
