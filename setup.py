from setuptools import setup

setup(name='ClusterPrepForDFT',
      version='0.1',
      description='Working on systematically generating clusters to do solid state modelling in a molecular DFT code, e.g. NWChem',
      author='WEspec',
      packages=['ClusterPrepForDFT'],
      requires=['pymatgen']
     )