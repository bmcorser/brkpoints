
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'version': '0.1',
    'description': 'Fragile Site Finder',
    'author': 'mc',
    'install_requires': ['nose','PyYAML'],
    'packages': ['brkpoints'],
    'scripts': [],
    'name': 'brkpoints',
    'package_data': {'tests': ['data/*'],},
    'entry_points': {
        'console_scripts': [
            'brkpoints-prepare = brkpoints.prepare:main',
            'brkpoints-find = brkpoints.find:main',
            ],
        },
}

setup(**config)
