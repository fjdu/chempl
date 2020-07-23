`Chempl` is an astrochemical code meant to be easy to play with in a python environment.

# Requirement

- Up-to-date compilers for Fortran, C++, and Cython
- Python 3
- For running the examples in the Jupyter notebooks, Anaconda3 is recommended

# Compile

- For compiling the python wrapper, run

    `make chempl`

- For compiling the command line executable, run

    `make`

  This will create an executable with the default name `re` (meaning "rate equation")

Depending on your system configuration, you _may_ need to edit the line

    library_dirs=["./"]

in the file "setup.py", which specifies the location of the gfortran library (libgfortran.a).

You may try to find the location of "libgfortran.a" using the following command

    dirname `gfortran --print-file-name libgfortran.a`


# Examples

- For the python wrapper

  Before detailed documentation becomes available, you may look at this [Jupyter notebook](https://github.com/fjdu/chempl/blob/master/Examples-2020-07-19.ipynb).

  Comments and suggestions are welcome!

- For the command line executable, you can run an example model by

    `./re paths_test.dat`

# Publication

- A paper describing `Chempl` can be downloaded from arxiv:
  [Chempl: a playable package for modeling interstellar chemistry](https://arxiv.org/abs/2007.11294)
  - Citation (BibTeX):
    ```@ARTICLE{2020arXiv200711294D,
           author = {{Du}, Fujun},
            title = "{Chempl: a playable package for modeling interstellar chemistry}",
          journal = {arXiv e-prints},
         keywords = {Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2020,
            month = jul,
              eid = {arXiv:2007.11294},
            pages = {arXiv:2007.11294},
    archivePrefix = {arXiv},
           eprint = {2007.11294},
     primaryClass = {astro-ph.SR},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200711294D},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    ```

