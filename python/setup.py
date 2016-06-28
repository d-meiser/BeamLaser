from setuptools import setup, find_packages
from codecs import open
from os import path
from distutils.extension import Extension
from Cython.Build import cythonize


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='beamlaser',

    version='0.0.1',

    description='Beam laser simultions',
    long_description=long_description,

    url='https://github.com/d-meiser/beamlaser',

    author='Dominic Meiser',
    author_email='dmeiser79@gmail.com',

    license='GPLv3',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GPLv3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    # What does your project relate to?
    keywords='laser',

    packages=find_packages(), 
    ext_modules = [Extension("beamlaser", ["beamlaser.pyx"],
       include_dirs=["../include", "../build-guile-mpi/", "../mpich/include"],
       extra_compile_args=["-std=c99"])],

    install_requires=['cython', 'numpy'],
)

