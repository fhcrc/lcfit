# Requirements

## virtualenv

``` shell
# create a virtualenv
$ virtualenv venv

# install packages
$ pip install -r requirements.txt
$ pip install --egg SCons==2.3.0
```

## R packages

``` shell
# create user library directory
mkdir venv/lib/R

# install packages into user library
$ R_LIBS_USER=$PWD/venv/lib/R R -e "install.packages('ape')"
```

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

# Activate environment

``` shell
$ source env.sh
```

# Running simulations

``` shell
$ scons .
```
