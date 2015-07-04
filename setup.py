#! /usr/bin/python2.7

from distutils.core import setup

setup(name='astools',
      description='astools = Atomic Simulation Tools. For use with CRYSTAL'+ 
                  'and CASTEP datafiles',
      author='Eric',
      url='https://github.com/mehese/astools',
      version='0.0.1',
      license='GPLv3',
      platforms=['linux'],
      long_description=open('README').read(),
      # Prerequisites
      install_requires=['numpy'],
      # Package info
      packages=['astools']
     )
