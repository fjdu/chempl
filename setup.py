from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
Options.docstrings = True

extension1 = Extension(
    name="chempl",
    sources=["chempl.pyx"],
    libraries=["chempl", "gfortran"],
    library_dirs=["./"],
    include_dirs=["./"],
    depends=['libchempl.a', 'setup.py', 'makefile'],
    extra_compile_args=["-std=c++11"],
)

extension2 = Extension(
    name="myconsts",
    sources=["myconsts.pyx"],
    libraries=["myconsts"],
    library_dirs=["./"],
    include_dirs=["./"],
    depends=['setup.py', 'makefile'],
    extra_compile_args=["-std=c++11"],
)

setup(
    name="chempl",
    version='0.1',
    description='A playable astrochemical code',
    author='Fujun Du',
    author_email='fjdu@pmo.ac.cn fujun.du@gmail.com',
    ext_modules=cythonize([extension1, extension2], language_level = "3",
    gdb_debug=True)
)
