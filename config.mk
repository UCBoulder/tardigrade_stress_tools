# Micromorphic filter, a variationally based filter for DNS
#
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the common configuration file for all of the included makefiles

# C++ compiler
CXX=/opt/moose/gcc-7.3.0/bin/g++

# Flags for the C++ compiler
CFLAGS=-std=gnu++11 -Wall -ansi -pedantic -O3 -I. -fmax-errors=5

# Location of the Eigen library
EIGEN=-I/projects/nathanm/usr/local/include/eigen-git-mirror

# Set the root directory
ROOTDIR=/projects/nathanm/constitutive_models

# Add the location of the error_tools to the include and library
ERRORSOURCE = $(ROOTDIR)/error_tools/src/cpp/error_tools.cpp
ERRORHEADER = $(ROOTDIR)/error_tools/src/cpp/error_tools.h
INC=-I$(ROOTDIR)/error_tools/src/cpp
LIB=-L$(ROOTDIR)/error_tools/src/cpp

# Add the location of the vector_tools to the include and library
VECTORSOURCE = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.cpp
VECTORHEADER = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.h
INC+=-I$(ROOTDIR)/vector_tools/src/cpp
LIB+=-L$(ROOTDIR)/vector_tools/src/cpp

# Add the location of the constitutive_tools to the include and library
CTSOURCE = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.cpp
CTHEADER = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.h
INC+=-I$(ROOTDIR)/constitutive_tools/src/cpp
LIB+=-L$(ROOTDIR)/constitutive_tools/src/cpp

# The python command
PYTHON=/apps/anaconda3/bin/python
