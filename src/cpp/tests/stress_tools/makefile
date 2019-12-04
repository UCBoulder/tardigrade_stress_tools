# Stress tools library
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   November 25, 2019
#
# This is the makefile for test_stress_tools.cpp it tests the 
# constitutive modeling tools.

include ../../../../config.mk

# List of executables
EXECUTABLE=test_stress_tools
TESTOBJECT=test_stress_tools.o
TESTSOURCE=test_stress_tools.cpp

OBJECT=stress_tools.o
SOURCE=../../stress_tools.cpp
HEADER=../../stress_tools.h

CTOBJECT=constitutive_tools.o
ERROROBJECT=error_tools.o

OBJECTS = $(TESTOBJECT)
OBJECTS += $(OBJECT) $(ERROROBJECT) $(CTOBJECT)

INC+=-I../..
LOCAL_CFLAGS=-ggdb

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(TESTOBJECT) $(OBJECT) $(ERROROBJECT) $(CTOBJECT)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11

$(TESTOBJECT): $(TESTSOURCE) $(HEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $< -c -o $@ $(CFLAGS) -lstdc++ -std=gnu++11

$(OBJECT): $(SOURCE) $(HEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

$(CTOBJECT): $(CTSOURCE) $(CTHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

#$(VECTOROBJECT): $(VECTORSOURCE) $(VECTORHEADER)
#	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

clean:
	rm -f $(EXECUTABLE)
	rm $(OBJECTS)

.PHONY: all clean
