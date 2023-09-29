.. _changelog:


#########
Changelog
#########

******************
0.5.4 (unreleased)
******************

******************
0.5.3 (2023-09-29)
******************

Internal Changes
================
- Add draft GitHub release action (:merge:`78`). By `Kyle Brindley`_.

******************
0.5.2 (2023-09-20)
******************
- Clean up conda package CI files after conda-build (:issue:`23`, :merge:`50`, :issue:`25`, :merge:`53`). 
  By `Sergio Cordova`_.

Internal Changes
================

******************
0.5.1 (2023-07-24)
******************

Breaking changes
================
- Change project, package, and namespace from 'solver tools' to 'tardigrade solver tools' (:issue:`15`, :merge:`48`). By
  `Kyle Brindley`_.

Internal Changes
================
- Clean up conda-build recipe (:issue:`22`, :merge:`46`). By `Kyle Brindley`_.
- Help CMake find the correct Python executable for conda-build on osx-arm64 (:merge:`47`). By `Kyle Brindley`_.
- Remove compiler as a runtime dependency. The OS-correct standard library package is added as a depedency by
  conda-build (:merge:`48`). By `Kyle Brindley`_.
- Build stdlib variants instead of compiler variants (:merge:`49`). By `Kyle Brindley`_.

******************
0.4.1 (2023-06-20)
******************

Breaking changes
================
- Deploy to the Conda environment preferred ``lib`` directory instead of the CMake linux default ``lib64`` (:issue:`21`,
  :merge:`44`). By `Kyle Brindley`_.

******************
0.3.1 (2023-04-03)
******************

Breaking Changes
================
- Require c++17 packages (:issue:`19`, :merge:`40`). By `Kyle Brindley`_.

Internal Changes
================
- Add a GCC 11 conda package variant (:issue:`16`, :merge:`33`). By `Kyle Brindley`_.
- Add the Sphinx target (:issue:`17`, :merge:`34`). By `Kyle Brindley`_.
- Force CI environment to build consistently from conda-forge (:merge: `35`). By `Nathan Miller`_.
- Prefer project-wide compiler options and remove ``-ansi`` to get consistent application of
  c++17 (:merge: `36`). By `Nathan Miller`_.
- Add a GCC 10 conda package variant (:issue:`18`, :merge:`37`). By `Sergio Cordova`_.
- Updates for parentheses/braces/brackets to match style guide (:merge:`38`). By `Kyle Brindley`_.
- Updated interface to the gradient of the determinant of a matrix w.r.t. the matrix (:merge:`43`). By `Nathan Miller`_.

******************
0.2.8 (2023-02-28)
******************

New Features
============
- Add an option to construct the stiffness tensor from the full 81 components (:issue:`11`, :merge:`23`). By `Kyle
  Brindley`_.
- Add an option to construct and rotate the stiffness tensor (:issue:`12`, :merge:`24`). By `Kyle Brindley`_.
- Add an energy and derivatives overload that accepts an Euler angle rotation for the stiffness matrix (:issue:`14`,
  :merge:`27`). By `Kyle Brindley`_.

Internal Changes
================
- Update minimum version requirements for ``tardigrade_vector_tools`` dependency (:merge:`25`). By `Kyle Brindley`_.
- Project configuration and conda build recipe changes to allow macOS builds and conda-build test stage (:merge:`16`).
  By `Kyle Brindley`_.
- Remove depreciated shell build script and documentation references in preference to direct cmake commands
  (:issue:`15`, :merge:`28`). By `Kyle Brindley`_.
- Update minimum ``tardigrade_vector_tools`` version requirement. By `Kyle Brindley`_.
- Fall back to aea-beta environment when project's CI environment doesn't exist (:merge:`32`). By `Kyle Brindley`_.

******************
0.2.7 (2022-12-21)
******************

New Features
============
- Added the computation of the gradients with respect to the previous parameter values for linear elasticity
  (:merge:`20`). By `Nathan Miller`_.

******************
0.2.6 (2022-12-16)
******************

New Features
============
- Add linear elasticity submodule from asp (:issue:`5`, :merge:`14`). By `Kyle Brindley`_.
- Add fully anisotropic, orthotropic, transverse isotropic, and cubic linear elasticity (:issue:`6`, :merge:`15`). By
  `Kyle Brindley`_.

Bug Fixes
=========
- Build and install a single shared library to help downstream projects find the full namespace (:issue:`10`,
  :merge:`17`). By `Kyle Brindley`_.

Internal Changes
================
- Updating framework to current cpp_stub standard (:merge:`12`). By `Nathan Miller`_.
- Removing additional errors preventing deploying the framework (:merge:`13`). By `Nathan Miller`_.
- Remove deprecated engineering constants stiffness tensor function interface (:issue:`9`, :merge:`16`). By `Kyle
  Brindley`_.

******************
0.2.5 (2022-03-21)
******************
