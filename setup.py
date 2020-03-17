from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
Options.docstrings = True

extension = Extension(
    name="chemplay",
    sources=["chemplay.pyx"],
    libraries=["chemplay", "gfortran"],
    library_dirs=["./", "/usr/local/Cellar/gcc/9.2.0_1/lib/gcc/9/"],
    include_dirs=["./"],
    depends=['libchemplay.a', 'setup.py', 'makefile'],
    extra_compile_args=["-std=c++11"],
)

setup(
    name="chemplay",
    version='0.0',
    description='A playable astrochemical code',
    author='Fujun Du',
    author_email='fjdu@pmo.ac.cn fujun.du@gmail.com',
    ext_modules=cythonize([extension], language_level = "3",
    gdb_debug=True)
)
