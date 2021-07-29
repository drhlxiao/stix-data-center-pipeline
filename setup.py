import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='stix',
    version='1.3',
    description='STIX python data parser',
    url='https://github.com/i4Ds/STIX-python-data-parser',
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'PyQt5', 'PyQtChart', 'scipy', 'pymongo', 'python-dateutil',
                      'xmltodict', 'spiceypy','qtconsole'],
    extras_require={
            'dev':  ["pytest", "pycodestyle", "pydocstyle",  "flake8"],
    },
    entry_points={
        'console_scripts': [
            'stix-parser= stix.apps.parser:main',
            'stix-quicklook= stix.apps.stix_quicklook:main',
        ],
        'gui_scripts': [
            'stix-parser-gui= stix.apps.stix:main',
        ]
    },
    python_requires='>=3.6'
)
