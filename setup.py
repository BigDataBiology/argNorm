from setuptools import setup
from setuptools import find_packages

# change this.
NAME = "argNorm"
AUTHOR = "Hui Chong"
EMAIL = "huichong.me@gmail.com"
URL = "https://github.com/BigDataBiology/argNorm"
LICENSE = "MIT"
DESCRIPTION = "Fast ARG normalization by mapping to the ARO ontology."


if __name__ == "__main__":
    setup(
        name=NAME,
        version="0.0.1rc4",
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        license=LICENSE,
        description=DESCRIPTION,
        packages=['argNorm', 'argNorm.data'],
        include_package_data=True,
        package_dir={'argNorm': 'argNorm' },
        package_data={'argNorm': ['data/*.tsv']},
        install_requires=open("./requirements.txt", "r").read().splitlines(),
        long_description=open("./README.md", "r").read(),
        long_description_content_type='text/markdown',
        entry_points={
            "console_scripts": [
                "argnorm=argNorm.cli:main"
            ]
        },
        zip_safe=False,
        classifiers=[
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3.8",
            "Development Status :: 4 - Beta",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Natural Language :: English"
        ],
        python_requires=">=3.8"
    )
