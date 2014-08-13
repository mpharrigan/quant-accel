from setuptools import setup, find_packages

setup(name='maccelerator',
      version='0.2',
      author='Matthew Harrigan',
      packages=find_packages(),
      scripts=['scripts/maccel.py'],
      zip_safe=False,
      package_data={'maccelerator':['reference/*.*']})
