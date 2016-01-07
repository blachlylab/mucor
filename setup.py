from setuptools import setup, find_packages

setup(name="mucor",
    version="1.5",
    description="Genomic Variant Aggregation and Mutation Correlation",
    author="James Blachly",
    author_email="james blachly at gmail com",
    license="GNU GPLv3",
    url="http://github.com/blachlylab/mucor",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "pandas", "HTSeq"],
    keywords = "sequencing, VCF, mutation, variant",
    classifiers = [ "Development Status :: 5 - Production/Stable",
                    "Environment :: Console",
                    "Intended Audience :: Healthcare Industry",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering :: Bio-Informatics" ]
    )
