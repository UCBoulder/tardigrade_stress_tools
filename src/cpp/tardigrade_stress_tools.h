/**
  ******************************************************************************
  * \file tardigrade_stress_tools.h
  ******************************************************************************
  * A collection of tools which implement and solve stress-strain relationships
  * in such a way to enable more rapid development of constitutive models which
  * have capabilities which may not be contained within a collection of
  * of constitutive models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_STRESS_TOOLS_H
#define TARDIGRADE_STRESS_TOOLS_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_error_tools.h>
#include<tardigrade_constitutive_tools.h>
#include<linear_elasticity.h>

namespace tardigradeStressTools{

    typedef tardigradeConstitutiveTools::floatType floatType; //!< Define the float values type.
    typedef tardigradeConstitutiveTools::floatVector floatVector; //!< Define a vector of floats
    typedef tardigradeConstitutiveTools::floatMatrix floatMatrix; //!< Define a matrix of floats

    void calculateMeanStress( const floatVector &stress, floatType &meanStress );

    floatType calculateMeanStress( const floatVector &stress );

    void calculateMeanStress( const floatMatrix &stress, floatType &meanStress );

    floatType calculateMeanStress( const floatMatrix &stress );

    void calculateMeanStress( const floatVector &stress, floatType &meanStress, floatVector &jacobian );

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric );

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric, floatMatrix &jacobian );

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric, floatMatrix &jacobian );

    floatVector calculateDeviatoricStress( const floatVector &stress );

    floatVector calculateDeviatoricStress( const floatVector &stress, floatMatrix &jacobian );

    void calculateVonMisesStress( const floatVector &stress, floatType &vonMises );

    floatType calculateVonMisesStress( const floatVector &stress );

    void calculateVonMisesStress( const floatVector &stress, floatType &vonMises, floatVector &jacobian );

    void druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B, floatType &dpYield );

    void druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam, floatType &dpYield );

    floatType druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B );

    floatType druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam );

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield );

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield );

    floatType druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B );

    floatType druckerPragerSurface( const floatVector &stress, const floatVector &dpParam );

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian );

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian );

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &djacobiandstress, const floatType tol = 1e-9 );

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress, const floatType tol = 1e-9 );

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &djacobiandstress );

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress );

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, const floatType tol = 1e-9 );

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, const floatType tol = 1e-9 );

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatVector &unitDirectionJacobian );

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian );

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatVector &unitDirectionJacobian );

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian );

//    void linearViscoelasticity( const floatType &currentTime, const floatType &currentStrain,
//                                   const floatType &previousTime, const floatType &previousStrain,
//                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
//                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
//                                   const floatType &alpha,
//                                   floatVector &stress, floatVector &currentStateVariables );
//
//    void linearViscoelasticity( const floatType &currentTime, const floatType &currentStrain,
//                                   const floatType &previousTime, const floatType &previousStrain,
//                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
//                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
//                                   const floatType &alpha,
//                                   floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain,
//                                   floatVector &dstressdrateModifier );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables,
                                   floatMatrix &dstressdstrain, floatVector &dstressdrateModifier );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables,
                                   floatMatrix &dstressdstrain, floatVector &dstressdrateModifier,
                                   floatMatrix &dstressdPreviousStrain, floatVector &dstressdPreviousRateModifier,
                                   floatMatrix &dstressdPreviousStateVariables,
                                   floatMatrix &dStateVariablesdStrain, floatVector &dStateVariablesdRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStrain, floatVector &dStateVariablesdPreviousRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStateVariables );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &dStress, floatVector &stress, floatVector &currentStateVariables );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &dStress, floatVector &stress, floatVector &currentStateVariables,
                                   floatMatrix &dstressdstrain, floatVector &dstressdrateModifier );

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &dStress, floatVector &stress, floatVector &currentStateVariables,
                                   floatMatrix &dstressdstrain, floatVector &dstressdrateModifier,
                                   floatMatrix &dstressdPreviousStrain, floatVector &dstressdPreviousRateModifier,
                                   floatMatrix &dstressdPreviousStateVariables,
                                   floatMatrix &dStateVariablesdStrain, floatVector &dStateVariablesdRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStrain, floatVector &dStateVariablesdPreviousRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStateVariables );

    void volumetricNeoHookean( const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress );

    void volumetricNeoHookean( const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ );

    void volumetricNeoHookean( const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress );

    void volumetricNeoHookean( const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ );

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatType n, floatType &p );

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p );

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatType n, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA );

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA );

    void linearHardening( const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value );

    void linearHardening( const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value, floatVector &valueJacobian );

    void computeJaumannStiffnessTensor( const floatVector &cauchyStress, const floatVector &currentDeformationGradient,
                                            const floatMatrix &dCauchydF, floatMatrix &C );
}

#ifdef TARDIGRADE_HEADER_ONLY
    #include "tardigrade_stress_tools.cpp"
#endif

#endif
