from setuptools import setup
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('pesci/version.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(name='pesci',
      version=main_ns['__version__'],
      description='Pretty Easy Single-cell Comparisons using ICC',
      url='http://github.com/eparey/pesci',
      author='Elise Parey',
      author_email='e.parey@ucl.ac.uk',
      license='TBD',
      packages=['pesci'],
      entry_points={'console_scripts': ['pesci = pesci.__main__:main']},
      install_requires=['coloredlogs', 'datatable', 'matplotlib', 'networkx', 'numpy>=2.0',
                        'pandas>=2.2.2', 'scanpy>=1.10.3', 'scipy', 'seaborn', 'tqdm'],
      zip_safe=False)