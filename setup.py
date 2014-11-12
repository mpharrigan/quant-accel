import os

from setuptools import setup, find_packages


REFDIR = 'maccelerator/reference'
if not os.path.exists(REFDIR):
    print("Making reference data")
    import make_reference_data

    make_reference_data.make_reference_data(REFDIR)

setup(name='maccelerator',
      version='0.2',
      author='Matthew Harrigan',
      packages=find_packages(),
      scripts=['scripts/maccel.py'],
      zip_safe=False,
      package_data={'maccelerator': ['reference/*.*']})
