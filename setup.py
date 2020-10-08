from setuptools import Extension
from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

extension = Extension("_atompack", sources=["extension/_atompack.c"])

setup(name="atompack",
      version="0.4.0",
      description="A flexible Python library for atomic structure generation.",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="Seaton Ullberg",
      author_email="seatonullberg@gmail.com",
      url="https://github.com/seatonullberg/atompack",
      license="MIT License",
      packages=find_packages(),
      ext_modules=[extension],
      extras_require={"dev": [
          "isort",
          "mypy",
          "pdoc3",
          "pyflakes",
          "pytest",
          "pytest-benchmark",
          "yapf",
      ]})
