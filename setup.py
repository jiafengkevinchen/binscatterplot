import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="binscatter",
    version="0.0.1",
    author="Jiafeng Chen",
    author_email="jchen@hbs.edu",
    description="Wrapper for R's binsreg in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jiafengkevinchen/binscatterplot",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
