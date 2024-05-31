from setuptools import setup, find_packages

setup(
    name = "AutoCoppersmith",
    version = "0.0.1",
    description = "AutoCoppersmith tools by Tearsjin",
    author = "Tearsjin",
    packages = find_packages(),
    install_requires = ["requests"],
    classifiers = [
        "Progamming Language :: Sagemath 9.5",
        "License:: OSI Approved :: MIT License"
    ]
)