import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="OctopusTools",
    version="0.0.1",
    author="xiao",
    author_email="xiao.chen@uni-jena.de",
    description="Some post data process of Octopus",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ChenXiao93/OctopusTools.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

