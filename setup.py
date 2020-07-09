from distutils.core import setup, Extension

import numpy as np

include_dirs = [
    "./extensions",
    np.get_include(),
]

extensions = [
    Extension(
        "_cell",
        sources=["./extensions/_cell.c"],
        include_dirs=include_dirs,
    ),
]

setup(
    name="atompack",
    version="0.1.0",
    description="A flexible Python library for atomic structure generation.",
    author="Seaton Ullberg",
    author_email="seatonullberg@gmail.com",
    url="https://github.com/seatonullberg/atompack",
    license="MIT License",
    ext_modules=extensions,
    packages=["atompack"],
)
