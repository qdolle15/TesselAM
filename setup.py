from setuptools import setup, find_packages

setup(
    name="TesselAM",
    version="1.0.1",
    description=(
        "TesselAM (Tessellation-based Additive Manufacturing) is a fast computational model "
        "for predicting competitive grain growth during directional solidification in metallic additive manufacturing. "
        "The name 'TesselAM' is inspired by the concept of tessellation, as the microstructure is represented "
        "through a set of seeds. Multiple seeds sharing the same identifier and orientation collectively define individual grains, "
        "enabling an efficient and physically consistent reconstruction of the evolving grain structure."
    ),
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Quentin Dollé",
    author_email="quentin.dolle@polytechnique.edu",
    url="https://github.com/qdolle15/TesselAM",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "tqdm",
        "matplotlib",
        "orix"
    ],
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    python_requires='>=3.8',
)