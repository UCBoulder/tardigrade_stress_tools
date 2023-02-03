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
#include<vector_tools.h>
#include<error_tools.h>
#include<constitutive_tools.h>

namespace stressTools{
namespace linearElasticity{

    typedef constitutiveTools::errorNode errorNode; //!< Redefinition for the error node
    typedef constitutiveTools::errorOut errorOut; //!< Redefinition for a pointer to the error node
    typedef constitutiveTools::floatType floatType; //!< Define the float values type.
    typedef constitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef constitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    errorOut formReferenceStiffnessTensor( const floatVector &parameters, floatMatrix &stiffnessTensor );

    errorOut formReferenceStiffnessTensor( const floatMatrix &directionCosines, const floatVector &parameters,
                                           floatMatrix &stiffnessTensor );

    errorOut formReferenceStiffnessTensor( const floatVector &bungeEulerAngles, const floatVector &parameters,
                                           floatMatrix &stiffnessTensor );

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy );

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress );

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi );

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                             floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 );

    errorOut evaluateEnergy( const floatVector &bungeEulerAngles,
                             const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                             floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 );

}  // linearElasticity
}  // stressTools

#endif
