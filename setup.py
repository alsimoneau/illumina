from setuptools import setup

setup(
    name='illum',
    version='2.1.21w46.1b-dev',
    py_modules=[
        'main',
        'alt_scenario_maker',
        'defineDomain',
        'extract-output-data',
        'find-failed-runs',
        'hdf_convert',
        'Illuminutils',
        'init_run',
        'make_inputs',
        'make_zones',
        'make_lamps',
        'makeBATCH',
        'MultiScaleData',
        'pytools',
        'street_orientation'
    ],
    install_requires=[
        'Click<8',
        'progressbar2',
        'pyproj',
        'pyyaml',
        'numpy',
        'h5py',
        'pillow',
        'matplotlib',
        'scipy',
        'astropy',
        'pandas',
        'geopandas',
        'fiona',
        'osmnx',
        'GitPython'
    ],
    entry_points='''
        [console_scripts]
        illum=main:illum
    ''',
)
