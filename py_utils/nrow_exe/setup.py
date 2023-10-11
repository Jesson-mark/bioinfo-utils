from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("nrow_exe.py"),
    zip_safe=False
)

