from setuptools import setup

setup(
    name='ligandome',
    version='0.1.0',
    description='Ligandome analysis for TCR-DSE paper',
    packages=['ligandome',
              'ligandome.utils',
              'ligandome.database_exports',
              'ligandome.third_party_tools']
)