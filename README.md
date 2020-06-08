# stress\_tools

Tools for computing stress-strain behaviors. The stress functions should be 
able to return the stress and the jacobian w.r.t. the strain metric of 
interest.

Note: In order to use the Intel compiler one must run the following command 
in a bash prompt:
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

This is the same command that the abaqus command issues. It may be that 
this command will change on different platforms.

---

---

## Dependencies

### Make

These tools have several dependencies that must be available in the same parent
directory as this repo. 

* eigen: https://gitlab.com/libeigen/eigen
* constitutive\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/constitutive_tools
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools

### CMake

The project is transitioning from Make to Cmake. For cmake builds, Eigen must be
"installed" following the ``eigen/INSTALL`` instructions. The Eigen dependence
is easiest to resolve if eigen is installed in the default install directory.
However, if you don't have admin privileges, you can also insall Eigen to your
home directory in ``~/include`` (or possibly in ``~/.local/include``, but this
is untested by this project).

#### Non-admin Eigen install for solver_tools
[Reference](https://unix.stackexchange.com/questions/36871/where-should-a-local-executable-be-placed)

```
# sstelmo
ssh -X sstelmo.lanl.gov
# source Intel compilers
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux
# Create personal include file directory
$ pwd
/home/$USER
$ mkdir include
# Move to repository directory
$ cd /preferred/path/to/repos
# Example
$ pwd
/projects/$USER/e13repos
# Clone eigen
$ git clone https://gitlab.com/libeigen/eigen.git
$ cd eigen
$ git checkout 3.3.7
# Build eigen
$ mkdir build
$ cd build
$ export CXX=$(command -v icpc)
$ cmake3 .. -DCMAKE_INSTALL_PREFIX=/home/$USER
$ make install
```

---

---

## Building the documentation

> **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
> API is much more useful

A build script has been created for convenience, ``new_build.sh``. It will build
everything including the library binary, the test binary, and the documentation.
This is the same build script used by ``jenkins_build.sh`` for CI builds and
testing.

### sstelmo

1) Activate the correct python environment

```
$ source /apps/anaconda/5.0.1-python-3.6/bin/activate
$ source activate /projects/python/release-cpp
```

2) Create the build directory and move there

```
$ pwd
/path/to/solver_tools/
$ mkdir build/
$ cd build/
```

3) Run cmake3 configuration

```
$ pwd
/path/to/solver_tools/build/
$ cmake3 ..
```

4) Build the docs

```
$ cmake3 --build docs
```

5) Documentation builds to: 

```
solver_tools/build/docs/sphinx/index.html
```

6) Display docs

```
$ pwd
/path/to/solver_tools/build/
firefox docs/sphinx/index.html &
```

7) While the Sphinx API is still a WIP, try the doxygen API

```
$ pwd
/path/to/solver_tools/build/
firefox docs/doxygen/html/index.html &
```
