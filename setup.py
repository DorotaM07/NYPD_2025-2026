from setuptools import setup, find_packages

setup(
    name="proteomika",
    version="1.0.0",
    description="Wstępna analiza proteomiki.",
    python_requires=">=3.9",
    packages=["proteomika"],
    install_requires=[
        "numpy>=1.24",
        "pandas>=2.0",
        "plotly>=5.0",
        "matplotlib>=3.7",
        "matplotlib-venn>=0.11",
    ],
    entry_points={
        "console_scripts": [
            "proteomika=proteomika.command:main",
        ],
    },
)
