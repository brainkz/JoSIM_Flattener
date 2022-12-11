from setuptools import setup
setup(
    name='netlist_flattener',
    version='0.0.1',
    entry_points={
        'console_scripts': [
            'myscript=myscript:run'
        ]
    }
)
