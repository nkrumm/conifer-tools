from setuptools import setup

setup(name='conifer-tools/',
      version='0.1',
      description='A collection of helpful tools for use with CoNIFER',
      url='http://github.com/nkrumm/conifer-tools',
      author='Niklas Krumm',
      author_email='nkrumm@gmail.com',
      license='',
      packages=['conifer-tools/'],
      install_requires=[
          'numpy',
          'tables',
          'pandas',
          'rpy2',
          'scipy',
          'drmaa'
      ],
      zip_safe=False)