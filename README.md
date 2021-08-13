# CHEOPSim
Source code for the CHEOPS simulator

<h3>Requirements</h3>

Installation of CHEOPSim requires the CHEOPS common software framework provided in [common_sw](https://github.com/davefutyan/common_sw) to be installed first.
See the [README.md](https://github.com/davefutyan/common_sw#readme) for that repository for instructions.

In addition to [common_sw](https://github.com/davefutyan/common_sw) and its dependencies, the CHEOPSim package requires the following:

* Gnu Scientific Library: http://www.gnu.org/software/gsl/
* C++ client API for PostgreSQL: http://pqxx.org/development/libpqxx/ (CHEOPS deployment uses version 4.0.1)

The following python packages are also required:
* numpy	1.18.1 or later
* scipy	1.1.0 or later

<h3>Installation</h3>

With the above dependencies installed, CHEOPSim can be installed as follows:

* cd CHEOPSim
* make install

<h3>Execution</h3>

CHEOPSim requires a set of reference files as input. A default set of reference files are provided in the tarball [reference_data.tar.gz](https://github.com/davefutyan/CHEOPSim/releases/download/v1.0/reference_data.tar.gz) provided with ther release of this package.

The files from this tarball need to be moved to ${CHEOPS_SW}/install

The configuration of the simulator is defined in CHEOPSim/simulator/conf/runCHEOPSim.xml

CHEOPSim can be run as follows:

* cd CHEOPSim/simulator
* runCHEOPSim

The runCHEOPSim command can be run from the CHEOPSim/simulator directory (as above) in which case, the configuration defined in CHEOPSim/simulator/conf/runCHEOPSim.xml will be used, or it can be executed from any directory containing a configuration file named runCHEOPSim.xml
