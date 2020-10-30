import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="xespresso",
    version="1.0.0",
    description="Quantum ESPRESSO Calculator for Atomic Simulation Environment (ASE).",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/superstar54/xespresso",
    author="Xing Wang",
    author_email="xingwang1991@gmail.com",
    license="GPL",
    classifiers=[
        "License :: GPL-3.0 License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    packages=["xespresso"],
    include_package_data=True,
    install_requires=["ase", "numpy", "scipy", "matplotlib"],
    entry_points={
        "console_scripts": [
            "realpython=reader.__main__:main",
        ]
    },
)
