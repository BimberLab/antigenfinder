import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="antigen_finder",
    version="1.0.0",
    author="Ben Bimber",
    author_email="benjamse@ohsu.edu",
    description="A tool for identifying neoantigens from variant data and identifying putative epitopes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BimberLab/antigen_finder",
    packages=setuptools.find_packages(),
    package_data={'antigenfinder': []},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[""],
    python_requires=">=3.10",
)
