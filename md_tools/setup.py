"""
Copyright (C) 2014,2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext as _build_ext


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


ext_module = Extension(
    'md_libs.bonds',
    ['md_libs/bonds.pyx'],
    )
ext_module_rdf = Extension(
    'md_libs._rdf',
    ['md_libs/_rdf.pyx'],
)

setup(
    name='lab-tools',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_module, ext_module_rdf],
    install_requires=['numpy', 'Cython'],
    test_suite='tests',
    version='0.0.1'
)
