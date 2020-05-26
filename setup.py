import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MAPT-SEpapoulis", # Replace with your own username
    version="0.9.0",
    author="Spiridon Evangelos Papoulis",
    author_email="spapouli@vols.utk.edu",
    description="A python package for designing a PNA clamp to block PCR amplification",
    license='GPLv3+',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SEpapoulis/MAPT",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
