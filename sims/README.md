# Requirements

## virtualenv

``` shell
# create a virtualenv
$ virtualenv venv

# activate virtualenv
$ source venv/bin/activate

# install packages
$ pip install -r requirements.txt
$ pip install --egg SCons==2.3.0
```

## R packages

``` shell
# restore packrat project
$ R --slave -e "0" --args --bootstrap-packrat
```

R packages are installed into the `packrat/lib` directory for use in the RStudio project.

## Bio++ libraries

installed Bio++ 2.2.0 via homebrew/science

## BppSuite

installed from [source](http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz) into the `venv` prefix, so it'll be in our `PATH` when the virtualenv is activated:

``` shell
$ wget http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz
$ tar xf bppsuite-2.2.0.tar.gz
$ cd bppsuite-2.2.0
$ cmake -DCMAKE_INSTALL_PREFIX=$PWD/../venv
$ make && make install
$ cd ..
```

# Running simulations

``` shell
$ scons .
```

The simulations can take several hours to complete, so you may wish to parallelize them by passing the `-j #` flag to `scons`, where `#` is the number of concurrent jobs to use.

# Running the analysis

The analysis isn't run as part of the SCons build above.
Instead, I've been running it in RStudio by sourcing the `bin/analyze_sims.R` file.
Before doing so, the working directory should be set to the `runs` directory.

TODO: run the analysis as part of the SCons build
