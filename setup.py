from setuptools import setup
from setuptools import find_packages

NAME = "argnorm"
AUTHOR = "See README"
EMAIL = "luispedro@big-data-biology.org"
URL = "https://github.com/BigDataBiology/argNorm"
LICENSE = "MIT"
DESCRIPTION = """
Normalize antibiotic resistance genes (ARGs) abundance tables (e.g., from metagenomics) by using the ARO ontology (developed by CARD).

"""


if __name__ == "__main__":
    setup(
        name=NAME,
        version="0.1.0",
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        license=LICENSE,
        description=DESCRIPTION,
        packages=['argnorm', 'argnorm.data'],
        include_package_data=True,
        package_dir={'argnorm': 'argnorm' },
        package_data={'argnorm': ['data/*.tsv']},
        install_requires=open("./requirements.txt", "r").read().splitlines(),
        long_description=open("./README.md", "r").read(),
        long_description_content_type='text/markdown',
        entry_points={
            "console_scripts": [
                "argnorm=argnorm.cli:main"
            ]
        },
        zip_safe=False,
        classifiers=[
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Development Status :: 4 - Beta",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Natural Language :: English"
        ],
    )
