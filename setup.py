from setuptools import setup, find_packages

setup(
    name="FAST-MMAM",
    version="0.1.0",
    description="FAST-MMAM: Fast simulation model for Metallic Microstructure Additive Manufacturing: \
        a computational model for competitive grain growth prediction during directional solidification \
        taking place during additive manufacturing process",
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
