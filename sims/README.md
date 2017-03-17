# Requirements

The following commands should all be run from within this directory (`sims`).


## virtualenv

``` shell
# create a virtualenv
virtualenv venv

# activate virtualenv
source venv/bin/activate

# install packages
pip install -r requirements.txt
pip install --egg SCons==2.3.0
```


## R packages

``` shell
# restore packrat project
R --no-restore --slave -e "0" --args --bootstrap-packrat
```

packrat installs the R packages into the `packrat/lib` directory.


## Bio++ and BppSuite

Running the simulations requires the `lcfit-compare` tool, so it should be built according to the instructions in the top-level `README.md`.
In addition, the tools `bppseqgen` and `bppml` from BppSuite 2.2.0 are required.

On Debian/Ubuntu:

``` shell
sudo apt-get install bppsuite
```

If BppSuite isn't available through your package manager, it can be installed from [source](http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz) into the `venv` prefix, so it'll be in your `PATH` when the virtualenv is activated:

``` shell
wget http://biopp.univ-montp2.fr/repos/sources/bppsuite/bppsuite-2.2.0.tar.gz
tar xf bppsuite-2.2.0.tar.gz
cd bppsuite-2.2.0
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../venv
make && make install
cd ..
```


# Running the simulations

``` shell
scons .
```

The simulations can take several hours to complete, so you may wish to parallelize them by passing the `-j #` flag to `scons`, where `#` is the number of concurrent jobs to use.


# Running the analysis

The simulation analysis isn't run as part of the SCons build above.
Instead, it can be run from the command line:

``` shell
R --no-restore --slave -e "setwd('runs'); source('../bin/analyze_sims.R')"
```

The analysis will produce the file `runs/measures.pdf` summarizing lcfit's performance.
