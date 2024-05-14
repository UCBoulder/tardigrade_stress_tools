/**
  ******************************************************************************
  * \file tardigrade_mass_change_deformation.h
  ******************************************************************************
  * The header file for the deformation that results from mass change.
  ******************************************************************************
  */

#ifndef TARDIGRADE_MASS_CHANGE_DEFORMATION_H
#define TARDIGRADE_MASS_CHANGE_DEFORMATION_H

#define USE_EIGEN

#include<vector>

namespace tardigradeStressTools{

    namespace massChangeDeformation{

        typedef double floatType;
        typedef std::vector< floatType > vector3d;
        typedef std::vector< floatType > secondOrderTensor;
        typedef std::vector< floatType > thirdOrderTensor;
        typedef std::vector< floatType > fourthOrderTensor;

    }

}

#include "tardigrade_mass_change_deformation.cpp"

#endif
