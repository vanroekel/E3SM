(omega-user-build)=

# CMake-based Omega Build

Omega's build system is constructed using the CMake build tool,
offering a strong foundation for managing the building process.

This page describes how to build Omega in general (outside Polaris or within
other workflows).

## If You Are Using Polaris

If you are already using Polaris to set up tasks or suites, you typically do
not need to build Omega manually. In that workflow, use Polaris auto-build
during `polaris setup` or `polaris suite`.

- Polaris quick start (developers):
  <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html>
- Building E3SM components:
  <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#building-e3sm-components>
- `omega_pr` suite:
  <https://docs.e3sm.org/polaris/main/users_guide/ocean/suites.html#omega-pr-suite>

Standalone/manual and Polaris-based build workflows are both supported.

## General Build Overview

The Omega build system has two modes: standalone and E3SM component.

For standalone build, the Omega build creates a generic E3SM case and
reads build configurations from the case. In E3SM component build, the
Omega build reads build configrations from the E3SM case created by
CIME build system. Following this, it configures CMake variables and
environment variables based on the computing system where the build is
taking place, as well as user input from the CMake command-line.

For the Omega build system to function, a Python interpreter is necessary.
The minimum version of CMake is 3.21.

## Standalone Build

CMake and Omega prefer an out-of-source build. This enables a user to build
and maintain multiple executables from the same source directory.
The standard practice is for a user to create a separate directory where
the build should take place and the commands below should be launched from
that directory.

To utilize a compiler name that's defined in the E3SM machine configuration
file (`config_machines.xml`), employ the `OMEGA_CIME_COMPILER` CMake variable,
as illustrated below. In some cases, you may also want to add `OMEGA_CIME_MACHINE`
to specify which system you intend to use. On systems where there is not a
default `PROJECT` defined, you can use `OMEGA_CIME_PROJECT` to specify an
account to use as a placeholder during the Omega build.
The values of `OMEGA_CIME_COMPILER` and `OMEGA_CIME_MACHINE` are defined in
"${E3SM}/cime\_config/machines/config\_machines.xml".
Omega requires some external libraries. Many of these are built automatically
from the E3SM distribution. However, the METIS, ParMETIS, and, optionally,
GKlib libraries must be built separately and the path must be supplied during
the cmake invocation as shown below. If a `OMEGA_METIS_ROOT` is not supplied,
it is assumed that both METIS and ParMETIS are installed in the same
`OMEGA_PARMETIS_ROOT` location. If a `OMEGA_GKLIB_ROOT` is not supplied, it is
assumed that both GKLIB and METIS are installed in the same `OMEGA_METIS_ROOT`
location. All three libraries should be static libraries.

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=nvidiagpu \
  -DOMEGA_CIME_MACHINE=pm-gpu \
  -DOMEGA_PARMETIS_ROOT=/path/to/parmetis \
  ${E3SM_HOME}/components/omega
```

Once the command completes successfully, several scripts will be created
in the build directory.

* omega\_env.sh   : load specific modules and set env. variables read from CIME
* omega\_build.sh : run `make` command after sourcing `omega_env.sh`
* omega\_run.sh   : run `./src/omega.exe` after sourcing `omega_env.sh`
* omega\_ctest.sh : run `ctest` after sourcing `omega_env.sh`

Run omega\_build.sh in the build directory to build Omega.

To employ a specific compiler, users can utilize the `OMEGA_CXX_COMPILER`
CMake variable, providing CMake with the compiler's path.

```sh
>> cmake \
  -DOMEGA_CXX_COMPILER=CC \
  ${E3SM_HOME}/components/omega
```

`OMEGA_CXX_COMPILER` overrides `OMEGA_CIME_COMPILER`.

To build Omega for GPU, add OMEGA\_ARCH to one of "CUDA" or "HIP".

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=nvidiagpu \
  -DOMEGA_CIME_MACHINE=pm-gpu \
  -DOMEGA_ARCH=CUDA \
  ${E3SM_HOME}/components/omega
```

To enable the ctest-based unittest option, include the `OMEGA_BUILD_TEST`
option as illustrated below.

```sh
>> cmake \
  -DOMEGA_BUILD_TEST=ON \
  ${E3SM_HOME}/components/omega
```

Once the `cmake` command succeeds, the directory where the command was
executed will contain the following files and directories.

```sh
>> ls -1 ${WORKDIR}
CMakeCache.txt
CMakeFiles
cmake_install.cmake
CTestTestfile.cmake
external
Makefile
omega_build.sh
omega_ctest.sh
omega_env.sh
omega_run.sh
src
test
```

To build the Omega, execute the `./omega_build.sh` command in the build directory.

If the build succeeds, the Omega library and executable are created in the
`src` sub-directory of the build directory.

To specify the output directory for saving the generated output files,
a user can include the `OMEGA_INSTALL_PREFIX` variable in the `cmake` command,
as shown below:

```sh
cmake \
  ... \
  -DOMEGA_INSTALL_PREFIX="/path/to/output/directory" \
  ${E3SM_HOME}/components/omega
  ...
```

Afterwards, a user can use the `make install` command to copy the Omega library
and executable to the `lib` and `bin` sub-directories under the directory
specified by `DOMEGA_INSTALL_PREFIX`.

The `./omega_ctest.sh` command runs Omega unit tests. To run the tests, MPI
parallel job launcher such as SLURM srun should be available. You may first
get allocation of an interactive computing node or use batch system.
In addition, you must either copy (or soft link) a mesh file into the same
directory as the unit test executables or specify in a configuration file
(not implemented yet) the path to that mesh input file. For detailed instructions
on how to obtain appropiate test meshes see {ref}`omega-dev-quick-start-getting-meshes`.

For example, if you are on an interactive computing node, you can run
Omega unit test by running `omega_ctest.sh` as shown below.

```sh
>> ./omega_ctest.sh
Test project <cmake working directory>
    Start 1: DATA_TYPES_TEST
1/1 Test #1: DATA_TYPES_TEST ..................   Passed    0.03 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.04 sec
```

(omega-user-cime-component-build)=

## E3SM Component Build

Omega can be built and run as the `ocn` component of a coupled E3SM simulation,
using a compset where Omega is the ocean component (see
`components/omega/cime_config/config_compsets.xml` for supported compsets).
Creating such a case and running `case.build` builds Omega
automatically: CIME's `components/omega/cime_config/buildlib_cmake` configures Omega with
`OMEGA_BUILD_MODE=E3SM` so that the coupler (c.f. standalone) driver is built.

(omega-cime-env)=

### Sourcing the CIME Environment

Omega uses YAML configuration files instead of traditional Fortran namelists.
When `case.setup` or `case.submit` runs, CIME executes Omega's
`components/omega/cime_config/buildnml` in-process within the Python interpreter
running CIME to generate `omega.yml`. Because `buildnml` depends on `pyyaml`,
which is not part of the Python standard library, `pyyaml` must be importable in
that Python environment before invoking any CIME commands.

Before running any CIME commands (`create_newcase`, `case.setup`, `case.build`,
`case.submit`, or `create_test`), source the CIME environment loader for your
machine:

| Machine | Environment Script |
| :--- | :--- |
| aurora | `source /lus/flare/projects/E3SMinput/soft/e3sm-unified/load_latest_cime_env.sh` |
| chrysalis | `source /lcrc/soft/climate/e3sm-unified/load_latest_cime_env.sh` |
| frontier | `source /ccs/proj/cli115/software/e3sm-unified/load_latest_cime_env.sh` |
| perlmutter (`pm-cpu`, `pm-gpu`) | `source /global/common/software/e3sm/anaconda_envs/load_latest_cime_env.sh` |

For a broader discussion of Python environment requirements and proposals for
a unified E3SM development and test environment, see
[E3SM Discussion #8906](https://github.com/E3SM-Project/E3SM/discussions/8906)
and the linked design document.

### Running a Coupled Case

Once the CIME environment is sourced, you can create, configure, build, and run
a coupled case from the root of your E3SM clone:

```bash
# 1. Source the CIME environment (e.g. on Chrysalis)
source /lcrc/soft/climate/e3sm-unified/load_latest_cime_env.sh

# 2. Navigate to CIME scripts and create a case
cd cime/scripts
./create_newcase --case /path/to/case --compset COMEGA-IAF --res T62_oQU240 --mach <machine>

# 3. Setup the case (generates omega.yml via buildnml)
cd /path/to/case
./case.setup

# 4. Build the case
./case.build

# 5. Submit the simulation
./case.submit
```

Supported compsets include `COMEGA-IAF` and `COMEGA-JRA1p5` (data atmosphere and
data sea ice with Omega ocean), as well as `GOMEGA-IAF` and `GOMEGA-JRA1p5`
(with active MPAS-Seaice). Model runtime options can be customized in
`user_nl_omega` prior to `case.setup` (see {ref}`omega-user-config`).
