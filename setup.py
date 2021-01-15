from setuptools import setup, find_packages

install_requires = ['numpy', 'astropy>=1.2, !=1.3.3', 'astroquery',
                    'scipy>=0.17.0', 'regions',
                    'tables', 'ipywidgets', 'astrowidgets', 'healpy']

extras = {'doc' : ['sphinx',  'sphinx-markdown-tables', 'sphinx-argparse',
                   'sphinx_rtd_theme']}

setup(
    # package description and version
    name="hetdex_api",
    version="0.8.3",
    author="The HETDEX Collaboration",
    author_email='erin@astro.as.utexas.edu',
    url='https://github.com/HETDEX/hetdex_api',
    download_url='https://github.com/HETDEX/hetdex_api/archive/0.8.2.tar.gz',
    description="Tools to deal with HETDEX data releases",

    # list of packages and data
    packages=find_packages(),
    include_package_data=True,
    # don't zip when installing
    zip_safe=False,

    # dependences
    install_requires=install_requires,
    extras_require=extras,

    classifiers=["Development Status :: 3 - Alpha",
                 "Environment :: Console",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License (GPL)",
                 "Operating System :: Unix",
                 "Programming Language :: Python :: 3.7",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 "Topic :: Utilities",
                 ],

    entry_points = {
                    "console_scripts" : [
                        'plot_completeness = hetdex_api.flux_limits.sensitivity_cube:plot_completeness',
                        'plot_completeness_versus_wl = hetdex_api.flux_limits.sensitivity_cube:plot_completeness_versus_wl',
                        'collapse_combine_sensitivity = hetdex_api.flux_limits.collapse_cubes:collapse_datacubes_command',
                        'biweight_fluxlims_hdf5 = hetdex_api.flux_limits.collapse_cubes:return_biwt_cmd',
                        'add_sensitivity_cube_to_hdf5 =  hetdex_api.flux_limits.hdf5_sensitivity_cubes:add_sensitivity_cube_to_hdf5',
                        'extract_sensitivity_cube = hetdex_api.flux_limits.hdf5_sensitivity_cubes:extract_sensitivity_cube',
                        'hetdex_get_spec = hetdex_tools.get_spec:main',
                        'hetdex_get_spec2D = hetdex_tools.get_spec2D:main',
                        'hetdex_get_shots = hetdex_tools.get_shots_of_interest:main'
                     ]
                   },

)
