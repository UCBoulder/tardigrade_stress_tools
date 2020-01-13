/*
===============================================================================
|                               stress_tools.h                                |
===============================================================================
| A collection of tools which implement and solve stress-strain relationships |
| in such a way to enable more rapid development of constitutive models which |
| have capabilities which may not be contained within a collection of         |
| of constitutive models.                                                     |
===============================================================================
*/

#ifndef STRESS_TOOLS_H
#define STRESS_TOOLS_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<constitutive_tools.h>

namespace stressTools{

    typedef constitutiveTools::errorNode errorNode; //!Redefinition for the error node
    typedef constitutiveTools::errorOut errorOut; //!Redefinition for a pointer to the error node
    typedef constitutiveTools::floatType floatType; //!Define the float values type.
    typedef constitutiveTools::floatVector floatVector; //! Define a vector of floats
    typedef constitutiveTools::floatMatrix floatMatrix; //!Define a matrix of floats

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress);

    floatType calculateMeanStress(const floatVector &stress);

    errorOut calculateMeanStress(const floatMatrix &stress, floatType &meanStress);

    floatType calculateMeanStress(const floatMatrix &stress);

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress, floatVector &jacobian);

    errorOut calculateDeviatoricStress(const floatVector &stress, floatVector &deviatoric);

    floatVector calculateDeviatoricStress(const floatVector &stress);

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises);

    floatType calculateVonMisesStress(const floatVector &stress);

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises, floatVector &jacobian);

    errorOut druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B, floatType &dpYield);

    errorOut druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam, floatType &dpYield);

    floatType druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B);

    floatType druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam);

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield);

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield);

    floatType druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B);

    floatType druckerPragerSurface(const floatVector &stress, const floatVector &dpParam);

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian);

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian);

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress);

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress);

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection);

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection);

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian);

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian);

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain, 
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier, 
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables);

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain, 
                                   const floatType &previousTime, const floatVector &previousStrain,               
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,                        
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain,
                                   floatVector &dstressdrateModifier);

    errorOut volumetricNeoHookean(const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress);

    errorOut volumetricNeoHookean(const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ);

    errorOut volumetricNeoHookean(const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress);

    errorOut volumetricNeoHookean(const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ);

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatType n, floatType &p);

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatType n, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA);

    errorOut linearHardening(const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift, 
                             floatType &value);
}

#endif
