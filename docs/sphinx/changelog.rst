.. _changelog:


#########
Changelog
#########

******************
0.2.5 (2022-12-21)
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
