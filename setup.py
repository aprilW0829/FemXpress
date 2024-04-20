import setuptools
from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FemXpress",
    version="v1.0.0",
    author="Xin Wang",
    author_email="wangxin970829@163.com",
    description="Tool that classifies the source of X chromosome inactivation of each female cell",
    long_description=long_description,
    url='https://github.com/wangxin970829/FemXpress',
    include_package_data=True,
    package_data={
    'FemXpress':['scripts/*']
    },
    #long_description_content_type="",
    #url="",
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    packages=setuptools find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL v3 or later",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)

