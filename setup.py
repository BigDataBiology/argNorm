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
        version="0.0.1rc1",
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        license=LICENSE,
        description=DESCRIPTION,
        packages=find_packages(),
        package_dir={'ArgNorm': 'ArgNorm'},
        include_package_data=True,
        install_requires=open("./requirements.txt", "r").read().splitlines(),
        long_description=open("./README.md", "r").read(),
        long_description_content_type='text/markdown',
        # change package_name to your package name.
        entry_points={
            "console_scripts": [
                "argnorm=ArgNorm.CLI:main"
            ]
        },
        # package_data={
        #     # change package_name to your package name.
        #     "data": ["./data/deeparg_ARO_mapping.tsv"],
        # },
        zip_safe=True,
        classifiers=[
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3.8",
            "Development Status :: 4 - Beta",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Natural Language :: English"
        ],
        python_requires=">=3.6"
    )
