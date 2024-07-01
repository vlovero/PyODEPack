import os
import subprocess
import sys
import numpy as np
from setuptools import Extension, setup
from setuptools.command.install import install as _install

BASE_DIR = os.path.dirname(__file__)
DEBUG = not True


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
                       extra_compile_args=["-std=c++17", "-march=native", "-Ofast" if not DEBUG else "-g"])

    setup(name="pyodepack",
          packages=["pyodepack"],
          version="0.1",
          author="Vincent Lovero",
          author_email="vllovero@ucdavis.edu",
          description="python package for integrating ODEs",
          install_requires=["numpy>=1.22",
                            "scipy>=1.10.1",
                            "setuptools>=61.0",
                            "wheel"],
          cmdclass={"install": install},
          include_package_data=True,
          package_data={"pyodepack": ["py.typed", "__init__.pyi"]},
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
