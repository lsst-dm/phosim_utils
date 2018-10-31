"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.phosim.utils


_g = globals()
_g.update(build_package_configs(
    project_name='phosim_utils',
    version=lsst.phosim.utils.version.__version__))
