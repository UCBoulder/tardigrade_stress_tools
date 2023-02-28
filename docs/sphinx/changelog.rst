.. _changelog:


#########
Changelog
#########

******************
0.2.9 (unreleased)
******************

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
- Update minimum version requirements for ``vector_tools`` dependency (:merge:`25`). By `Kyle Brindley`_.
- Project configuration and conda build recipe changes to allow macOS builds and conda-build test stage (:merge:`16`).
  By `Kyle Brindley`_.
- Remove depreciated shell build script and documentation references in preference to direct cmake commands
  (:issue:`15`, :merge:`28`). By `Kyle Brindley`_.

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
