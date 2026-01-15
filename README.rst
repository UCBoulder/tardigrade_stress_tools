.. targets-start-do-not-remove

.. _Doxygen: https://www.doxygen.nl/manual/docblocks.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _PEP-8: https://www.python.org/dev/peps/pep-0008/
.. _`gersemi`: https://github.com/BlankSpruce/gersemi
.. _`clang-tidy`: https://clang.llvm.org/extra/clang-tidy/
.. _`clang-format`: https://clang.llvm.org/docs/ClangFormat.html

.. targets-end-do-not-remove

#############
stress\_tools
#############

*******************
Project Description
*******************

.. project-brief-start-do-not-remove

Tools for computing stress-strain behaviors. The stress functions should be
able to return the stress and the jacobian w.r.t. the strain metric of
interest.

.. project-brief-end-do-not-remove

Information
===========

TODO

Developers
==========

* Nathan Miller Nathan.A.Miller@colorado.edu
* Kyle Brindley kbrindley@lanl.gov

************
Dependencies
************

.. dependencies-start-do-not-remove

Compilers
=========

* c++11 compiler (listed version number has been tested at some point)

  * g++ >= GNU 4.8.5

Executables
===========

* [CMake](https://cmake.org/cmake/help/v3.14/) >= 3.14
* [Doxygen](https://www.doxygen.nl/manual/docblocks.html) >= 1.8.5
* [LaTeX](https://www.latex-project.org/help/documentation/) >= 2017

Python Modules (for documentation)
==================================

For convenience, the minimal Python environment requirements for the documentation build are included in
``environment.txt``. A minimal anaconda environment for building the documentation can be created from an existing
anaconda installation with the following commands.

.. code-block:: bash

   $ conda create --name tardigrade_stress_tools-env --file environment.txt --channel file:///projects/aea_compute/aea-conda --channel conda-forge

You can learn more about Anaconda Python environment creation and management in the [Anaconda
Documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

C++ Libraries
=============

.. note::

   **NOTE: Non-admin installations for Eigen and Boost are no longer required.** This project is built and deployed
   against C++ libraries managed in Conda. See the Conda environment file and README discussion for non-admin environment
   management.

* [Eigen](https://eigen.tuxfamily.org/dox/) >= 3.3.7
* [BOOST](https://www.boost.org/doc/libs/1_53_0/) >= 1.53.0
* error\_tools: https://github.com/UCBoulder/tardigrade_error_tools
* vector\_tools: https://github.com/UCBoulder/tardigrade_vector_tools
* constitutive\_tools: https://github.com/UCBoulder/tardigrade_constitutive_tools

If not found on the current system or active Conda environment, all of the
``*_tools`` libraries are pulled from their git repos by branch name and built
with their respective cmake files as part of the cmake build for this project.

.. dependencies-end-do-not-remove

**************
Build and Test
**************

.. build-start-do-not-remove

This project is built with [CMake](https://cmake.org/cmake/help/v3.14/) and uses
[Sphinx](https://www.sphinx-doc.org/en/master/) to build the documentation with
[Doxygen](https://www.doxygen.nl/manual/docblocks.html) +
[Breathe](https://breathe.readthedocs.io/en/latest/) for the c++ API.

.. warning::

   **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
   API is much more useful

sstelmo
=======

1) Activate the shared development environment

   .. code-block:: bash

      $ module use /projects/aea_compute/modulefiles
      $ module load tardigrade_stress_tools-env

2) Build everything

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/
      $ mkdir build
      $ cd build
      $ cmake ..
      $ cmake --build . --target all

3) View test results

   .. code-block:: bash

      cat build/src/cpp/tests/results.tex

4) Display docs

   .. code-block:: bash

      # Sphinx
      firefox build/docs/sphinx/html/index.html &

      # Doxygen
      firefox build/docs/doxygen/html/index.html &

Local development
=================

In some cases it is not convenient to pull down every repository required but it may be desired that local
versions of the repository are used. An example of when this may be needed is if development is across
multiple libraries and is proceeding faster than collaborators can check in results. In this case, and
outside of developers no-one should need to do this, a version of the code using local repositories can be
built.

To perform in-source builds of upstream libraries, the active Conda environment can NOT include installed versions of
the upstream libraries to be built in-source with the current project. It is possible to mix sources with some upstream
libraries coming from the active Conda environment and others built in-source from a Git repository. Developers may
build minimal working Conda environments from the Python Modules discussion.

1) Build and activate a minimal Conda development environment

   .. code-block:: bash

       $ conda create --name tardigrade_stress_tools-env --file environment.txt --channel file:///projects/aea_compute/aea-conda --channel conda-forge
       $ conda activate tardigrade_stress_tools-env

2) Define convenience environment variables

   .. code-block:: bash

       $ tardigrade_error_tools=/path/to/my/tardigrade_error_tools
       $ tardigrade_error_tools_version=origin/dev
       $ tardigrade_vector_tools=/path/to/my/tardigrade_vector_tools
       $ tardigrade_vector_tools_version=origin/dev

3) Perform the initial configuration. Note that the environment variables are mutually independent. Each variable can be
   used alone or in arbitrary combinations. The default values are found in the root ``CMakeLists.txt`` file. The ``PATH``
   variables can accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_REPOSITORY`` option can accept.
   The ``GITTAG`` variables will accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_TAG`` option can accept.

   .. code-block:: bash

      # View the defaults
      $ grep _TOOLS_ CMakeLists.txt
      set(TARDIGRADE_ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of tardigrade_error_tools")
      set(TARDIGRADE_ERROR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of tardigrade_error_tools")
      set(TARDIGRADE_VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of tardigrade_vector_tools")
      set(TARDIGRADE_VECTOR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of tardigrade_vector_tools")

      $ Build against local directory paths and possible custom branch
      $ pwd
      /path/to/tardigrade_stress_tools
      $ mkdir build
      $ cd build
      $ cmake .. -DFETCH_SOURCE=LOCAL -DTARDIGRADE_ERROR_TOOLS_PATH=${tardigrade_error_tools} -DTARDIGRADE_VECTOR_TOOLS_PATH=${tardigrade_vector_tools}

4) Building the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build
      $ make


Building the documentation
==========================

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/
      $ mkdir build/
      $ cd build/

3) Run cmake3 configuration

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build/
      $ cmake3 ..

4) Build the docs

   .. code-block:: bash

      $ cmake3 --build docs

5) Documentation builds to:

   .. code-block:: bash

      tardigrade_stress_tools/build/docs/sphinx/index.html

6) Display docs

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build/
      $ firefox docs/sphinx/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build/
      $ firefox docs/doxygen/html/index.html &

.. build-end-do-not-remove

*******************
Install the library
*******************

Build the entire before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build
      $ cmake3 --build .

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_stress_tools/build
      $ cmake --install . --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install . --prefix /home/$USER/.local

      # Example install to conda environment
      $ conda activate my_env
      $ cmake --install . --prefix ${CONDA_DEFAULT_ENV}

***********************
Contribution Guidelines
***********************

.. contribution-start-do-not-remove

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

`Sphinx`_ reads in docstrings and other special portions of the code as
reStructured text. Developers should follow styles in this `Sphinx style guide
<https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#>`_.

Style Guide
===========

This project uses the `gersemi`_ CMake linter. The CI style guide check runs the following command

.. code-block:

   $ gersemi CMakeLists.txt src/ docs/ --check

and any automatic fixes may be reviewed and then applied by developers with the following commands

.. code-block:

   $ gersemi CMakeLists.txt src/ docs/ --diff
   $ gersemi CMakeLists.txt src/ docs/ --in-place

This project enforces its style using `clang-tidy`_ and `clang-format`_ as configured with the
`.clang-format` and `.clang-tidy` files in the root directory. The formatting of the project can be
checked using `clang-tidy`_ by first configuring the project using

.. code-block:

   $ cmake -S . -B build ... -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

where `...` are the other configuration flags specified. After this clang-tidy can be run on the
full project from the source directory via

.. CAUTION::
    Commit all changes prior to running the clang tidy command. This will edit all source files.

.. code-block:

   $ run-clang-tidy -config-file=.clang-tidy -p build -extra-arg="-mno-sse2"

The formatting can be checked using `clang-format`_ by running

.. code-block:

   $ cmake -S . -B build ...
   $ cmake --build build --target cpp-format-check

which will indicate if the formatting is correct. The c++ files can be re-formatted to match the
style guidance by running

.. CAUTION::
    Commit all changes prior to running the format command. This will edit all source files.

.. code-block

   $ cmake --build build --target cpp-format

If the style is not constrained by the above, it should be inferred by the surrounding code.
Wherever a style can't be inferred from surrounding code this project falls back to `PEP-8`_-like
styles the exceptions to the notional PEP-8 fall back:

1. `Doxygen`_ style docstrings are required for automated, API from source documentation.

.. contribution-end-do-not-remove
