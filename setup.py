from setuptools import setup, Extension

import numpy as np

with open("README.md", "r") as f:
    long_description = f.read()

include_dirs = [
    "./extensions/",
    np.get_include(),
]

extensions = [
    Extension(
        "_cell",
        sources=["./extensions/cell.c"],
        include_dirs=include_dirs,
    ),
    Extension(
        "_pbc",
        sources=["./extensions/pbc.c"],
        include_dirs=include_dirs,
    ),
]

setup(name="atompack",
      version="0.3.0",
      description="A flexible Python library for atomic structure generation.",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="Seaton Ullberg",
      author_email="seatonullberg@gmail.com",
      url="https://github.com/seatonullberg/atompack",
      license="MIT License",
      ext_modules=extensions,
      packages=["atompack"],
      install_requires=["numpy", "scipy"],
      extras_require={"dev": [
          "pytest",
          "pytest-benchmark",
          "pdoc3",
          "isort",
          "yapf",
          "mypy",
          "pyflakes",
      ]})
