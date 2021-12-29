from setuptools import setup

setup(name='dnacycp',
      packages=['dnacycp'],
      version='0.0.1dev1',
      install_requires=[
      'numpy',
      'pandas',
      'tensorflow',
      'keras',
      'bio',
      'docopt'
      ],
      entry_points={
            'console_scripts': ['dnacycp-cli=dnacycp.cli:main']
      }
      )
