import os
import sys
import subprocess
import numpy as np
from setuptools import setup, Extension
from setuptools.command.install import install as _install


BASE_DIR = os.path.dirname(__file__)


class install(_install):
    def run(self):
        cmd = [sys.executable, __file__, "build"]
        process = subprocess.Popen(cmd)
        process.wait()
        _install.run(self)


if __name__ == "__main__":
    module = Extension("pyodepack._pyodepack",
                       sources=["pyodepack/src/_pyodepack.cpp"],
                       include_dirs=["pyodepack/include", np.get_include()],
                       extra_compile_args=["-std=c++17", "-march=native", "-Ofast"])

    setup(name="pyodepack",
          packages=["pyodepack"],
          version="0.1",
          author="Vincent Lovero",
          author_email="vllovero@ucdavis.edu",
          description="python package for integrating ODEs",
          install_requires=["setuptools", "numpy>=1.22", "scipy"],
          cmdclass={'install': install},
          include_package_data=True,
          package_data={'pyodepack': ['py.typed', '__init__.pyi']},
          classifiers=["Development Status :: 4 - Beta",
                       "Operating System :: MacOS",
                       "License :: OSI Approved :: MIT License",
                       "Programming Language :: C++",
                       "Programming Language :: Python :: 3.6",
                       "Programming Language :: Python :: 3.7",
                       "Programming Language :: Python :: 3.8",
                       "Programming Language :: Python :: 3.9",
                       "Programming Language :: Python :: 3.10"],
          ext_modules=[module])
