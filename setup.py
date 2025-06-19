from setuptools import setup, find_packages

setup(
    name="TesselAM",
    version="0.1.0",
    description="TesselAM (Tessellation-based Additive Manufacturing) is a fast computational model " + \
                "for predicting competitive grain growth during directional solidification in metallic additive manufacturing. " + \
                "The name 'TesselAM' is inspired by the concept of tessellation, as the microstructure is represented " + \
                "through a set of seeds. Multiple seeds sharing the same identifier and orientation collectively define individual grains, " + \
                "enabling an efficient and physically consistent reconstruction of the evolving grain structure.",
    author="Quentin Dollé",
    author_email="quentin.dolle@polytechnique.edu",
    url="https://",
    packages=find_packages(),  # Check an __init__.py file is in every sub directories
    install_requires=[
        "numpy",
        "scipy",
        "tqdm",
        "matplotlib"
    ],
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    python_requires='>=3.7',
)
