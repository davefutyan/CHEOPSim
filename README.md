# CHEOPSim
Source code for the simulator of the [CHEOPS](https://cheops.unibe.ch) space telescope mission. CHEOPSim is described in the following publication:
https://arxiv.org/abs/2001.05587

<h3>Requirements</h3>

Installation of CHEOPSim requires the CHEOPS common software framework provided in [common_sw](https://github.com/davefutyan/common_sw) to be installed first.
See the [README.md](https://github.com/davefutyan/common_sw#readme) for that repository for instructions.

In addition to [common_sw](https://github.com/davefutyan/common_sw) and its dependencies, the CHEOPSim package requires the following:

* Gnu Scientific Library (GSL): http://www.gnu.org/software/gsl/
* C++ client API for PostgreSQL: http://pqxx.org/development/libpqxx/ (CHEOPS deployment uses version 4.0.1)
* Eigen C++ template library for linear algebra: https://eigen.tuxfamily.org/ (CHEOPS deployment uses version 3.3.4)


The following python packages are also required:
* numpy	1.18.1 or later
* scipy	1.1.0 or later

<h3>Environment</h3>

The environment first needs to be set up as described in the
[README.md for common_sw](https://github.com/davefutyan/common_sw#readme).

In addition, CHEOPSim requires the following environment variable to be set:

    export CHEOPS_TESTDATA=${CHEOPS_SW}

The following updates are also required:
* Add include paths for GSL, pqxx and Eigen to the EXT_INC_DIRS environment variable
* Add library paths for GSL and pqxx to the EXT_LIB_DIRS and LD_LIBRARY_PATH environment variables

<h3>Installation</h3>
CHEOPSim can now be installed as follows:

    cd CHEOPSim
    make install

<h3>Execution</h3>

CHEOPSim requires a set of reference files as input. A default set of reference files are provided in the tarball [reference_data.tar.gz](https://github.com/davefutyan/CHEOPSim/releases/download/v1.0/reference_data.tar.gz) provided with the release of this package.

The contents of this tarball need to be moved to ${CHEOPS_SW}/install/resources

The configuration of the simulator is defined in CHEOPSim/simulator/conf/runCHEOPSim.xml

The parameters of the configuration are documented in the [CHEOPSim user guide](https://github.com/davefutyan/CHEOPSim/releases/download/v1.0/CHEOPSim_UserManual.pdf) provided with the release of this package.

The configuration file can be generated using the web interface provided in the [CHEOPSim_web](https://github.com/davefutyan/CHEOPSim_web) package.

CHEOPSim can be run as follows:

    cd CHEOPSim/simulator
    runCHEOPSim

The runCHEOPSim command can be run from the CHEOPSim/simulator directory (as above) in which case, the configuration defined in CHEOPSim/simulator/conf/runCHEOPSim.xml will be used, or it can be executed from any directory containing a configuration file named runCHEOPSim.xml
