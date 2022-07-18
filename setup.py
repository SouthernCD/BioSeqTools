# coding utf8
import setuptools
from seqparser.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="SeqParser",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A small tool for analyzing biological sequences",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/seqparser",

    entry_points={
        "console_scripts": ["SeqParser = seqparser.cli:main"]
    },    

    packages=setuptools.find_packages(),

    install_requires=[
        "toolbiox>=0.0.4",
        "bcbio-gff>=0.6.6",
        "biopython>=1.76",
        "pyfaidx>=0.5.5.2",
    ],

    python_requires='>=3.5',
)
