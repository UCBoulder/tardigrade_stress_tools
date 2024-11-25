/**
  ******************************************************************************
  * \file linear_elasticity.h
  ******************************************************************************
  * The header file for an implementation of quadratic energy form linear
  * elasticity. This serves as a test for the micro particle
  ******************************************************************************
  */

#ifndef LINEARELASTICITY_H
#define LINEARELASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_error_tools.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeStressTools{
namespace linearElasticity{

    typedef tardigradeConstitutiveTools::floatType floatType; //!< Define the float values type.
    typedef tardigradeConstitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef tardigradeConstitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    void formReferenceStiffnessTensor( const floatVector &parameters, floatMatrix &stiffnessTensor );

    void formReferenceStiffnessTensor( const floatMatrix &directionCosines, const floatVector &parameters,
                                           floatMatrix &stiffnessTensor );

    void formReferenceStiffnessTensor( const floatVector &bungeEulerAngles, const floatVector &parameters,
                                           floatMatrix &stiffnessTensor );

    void evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy );

    void evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress );

    void evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                         floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi );

    void evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                         floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                         floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 );

    void evaluateEnergy( const floatVector &bungeEulerAngles,
                         const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                         floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                         floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 );

}  // linearElasticity
}  // tardigradeStressTools

#endif
