/**
  * \file test_stress_tools.cpp
  *
  * Tests for stress_tools
  */

#include<stress_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_stress_tools
#include <boost/test/included/unit_test.hpp>

typedef constitutiveTools::errorOut errorOut;
typedef constitutiveTools::floatType floatType;
typedef constitutiveTools::floatVector floatVector;
typedef constitutiveTools::floatMatrix floatMatrix;

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer )
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( testCalculateMeanStress ){
    /*!
     * Test the mean stress calculation
     */

    //Initialize test values
    floatVector stressVector = { 1., 0., 0.,
                                 0., 1., 0.,
                                 0., 0., 1. };
    floatMatrix stressMatrix = { { 1., 0., 0. },
                                 { 0., 1., 0. },
                                 { 0., 0., 1. } };
    floatType meanStressExpected = 1.;
    floatVector jacobianVectorExpected = { 1./3.,       0.,       0.,
                                              0.,    1./3.,       0.,
                                              0.,       0.,    1./3. };

    //Initialize test output
    floatType meanStress;
    floatVector jacobianVector( stressVector.size( ) );

    //Test for correct mean stress calculation from stressVector with pointer output
    meanStress = 0.;
    errorOut result = stressTools::calculateMeanStress( stressVector, meanStress );
    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, meanStressExpected ) );

    //Test for correct mean stress calculation from stressVector without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress( stressVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, meanStressExpected ) );

    //Test for correct mean stress calculation from stressMatrix with pointer output
    meanStress = 0.;
    result = stressTools::calculateMeanStress( stressMatrix, meanStress );
    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, meanStressExpected ) );

    //Test for correct mean stress calculation from stressMatrix without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress( stressMatrix );
    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, meanStressExpected ) );

    //Test for correct mean stress and jacobian calculation
    meanStress = 0.;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    result = stressTools::calculateMeanStress( stressVector, meanStress, jacobianVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, meanStressExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

}

BOOST_AUTO_TEST_CASE( testCalculateDeviatoricStress ){
    /*!
     * Test the deviatoric stress calculation
     */
    floatVector stressVector = { 1., 0., 0.,
                                 0., 1., 0.,
                                 0., 0., 1. };
    floatVector expectedVector = { 0., 0., 0.,
                                   0., 0., 0.,
                                   0., 0., 0. };
    floatVector deviatoricVector( stressVector.size( ), 0. ), deviatoricVectorJ;
    floatMatrix jacobian, jacobian2;
    floatType eps = 1e-6;
    errorOut result;

    //Test computation of deviatoric tensor in row major format
    std::fill( deviatoricVector.begin( ), deviatoricVector.end( ), 0. );
    result = stressTools::calculateDeviatoricStress( stressVector, deviatoricVector );

    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals( expectedVector, deviatoricVector ) );

    std::fill( deviatoricVector.begin( ), deviatoricVector.end( ), 0. );

    deviatoricVector = stressTools::calculateDeviatoricStress( stressVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( expectedVector, deviatoricVector ) );

    //Test the computation of the jacobian
    result = stressTools::calculateDeviatoricStress( stressVector, deviatoricVectorJ, jacobian );

    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals( expectedVector, deviatoricVectorJ ) );

    //Test the jacobian
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps * fabs( stressVector[ i ] ) + eps;

        result = stressTools::calculateDeviatoricStress( stressVector + delta, deviatoricVectorJ, jacobian );

        BOOST_CHECK( ! result );

        floatVector grad = ( deviatoricVectorJ - deviatoricVector ) / delta[ i ];

        for ( unsigned int j=0; j<grad.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( grad[ j ], jacobian[ j ][ i ] ) );
        }

    }

    deviatoricVector = stressTools::calculateDeviatoricStress( stressVector, jacobian2 );

    BOOST_CHECK( vectorTools::fuzzyEquals( deviatoricVector, expectedVector ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( jacobian, jacobian2 ) );

}

BOOST_AUTO_TEST_CASE( testCalculateVonMisesStress ){
    /*!
     * Test the von Mises stress calculation
     */

    //Initialize test values
    floatVector stressVector = { 1., 1., 1.,
                                 1., 1., 1.,
                                 1., 1., 1. };
    floatVector jacobianVectorExpected = {    0., 1./2., 1./2.,
                                           1./2.,    0., 1./2.,
                                           1./2., 1./2.,    0. };
    floatType vonMisesExpected = 3.0;

    //Initialize test output
    floatType vonMises;
    floatVector jacobianVector( stressVector.size( ) );

    //Test computation of vonMises stress from row major stress tensor
    vonMises = 0.;
    stressTools::calculateVonMisesStress( stressVector, vonMises );
    BOOST_CHECK( vectorTools::fuzzyEquals( vonMises, vonMisesExpected ) );

    vonMises = 0.;
    vonMises = stressTools::calculateVonMisesStress( stressVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( vonMises, vonMisesExpected ) );

    //Test computation of vonMises stress and jacobian
    vonMises = 0.;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    stressTools::calculateVonMisesStress( stressVector, vonMises, jacobianVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( vonMises, vonMisesExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

}

BOOST_AUTO_TEST_CASE( testDruckerPragerSurface ){
    /*!
     * Test the Drucker-Prager yield criterion calculation.
     */

    //Declare test input variables
    floatVector stressVector = { 1., 1., 1.,
                                 1., 1., 1.,
                                 1., 1., 1. };
    floatType vonMises = 3.;
    floatType meanStress = 1.;
    floatType A = 1.;
    floatType B = 0.;
    floatVector dpParam = { A, B };

    floatType dpYieldExpected = 4.;
    floatVector jacobianVectorExpected = {  1./3.,  1./2.,  1./2.,
                                            1./2.,  1./3.,  1./2.,
                                            1./2.,  1./2.,  1./3. };
    floatVector unitDirectionVectorExpected = {  1./3.,  1./2.,  1./2.,
                                                 1./2.,  1./3.,  1./2.,
                                                 1./2.,  1./2.,  1./3. };
    unitDirectionVectorExpected /= sqrt( 1.5 + 1./3 );

    //Declare internal testing variables
    errorOut error;
    floatType eps;
    floatVector delta( stressVector.size( ), 0. );
    floatVector gradCol( stressVector.size( ), 0 );

    //Declare test output variables
    floatType dpYield;
    floatVector jacobianVector( stressVector.size( ), 0. );
    floatVector unitDirectionVector( stressVector.size( ), 0. );
    floatVector jacobianVectorJ( stressVector.size( ), 0. );
    floatMatrix djacobiandstress;
    floatMatrix unitDirectionJacobian;
    floatVector unitDirectionVectorJ( stressVector.size( ), 0. );

    //Test computation of DP yield criterion from vonMises and meanStress
    dpYield = 0;
    stressTools::druckerPragerSurface( vonMises, meanStress, A, B, dpYield );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    stressTools::druckerPragerSurface( vonMises, meanStress, dpParam, dpYield );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface( vonMises, meanStress, A, B );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface( vonMises, meanStress, dpParam );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    //Test computation of DP yield criterion from row major stress tensor
    dpYield = 0;
    stressTools::druckerPragerSurface( stressVector, A, B, dpYield );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    stressTools::druckerPragerSurface( stressVector, dpParam, dpYield );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface( stressVector, A, B );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface( stressVector, dpParam );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) );

    //Test computation of DP yield and jacobian from row major stress tensor
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    stressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    stressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

    //Test computation of DP yield, jacobian, and unit normal from row major stress tensor
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( unitDirectionVector.begin( ), unitDirectionVector.end( ), 0. );
    stressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, unitDirectionVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) &&
                 vectorTools::fuzzyEquals( unitDirectionVector, unitDirectionVectorExpected ) );

    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( unitDirectionVector.begin( ), unitDirectionVector.end( ), 0. );
    stressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector );
    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) &&
                 vectorTools::fuzzyEquals( unitDirectionVector, unitDirectionVectorExpected ) );

    //Test the computation of the DP yield, jacobian, and the derivative of the jacobian
    //w.r.t. the stress
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( jacobianVectorJ.begin( ), jacobianVectorJ.end( ), 0. );
    for ( unsigned int i=0; i<djacobiandstress.size( ); i++ ){
        std::fill( djacobiandstress[ i ].begin( ), djacobiandstress[ i ].end( ), 0. );
    }

    error = stressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, djacobiandstress );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

    eps = 1e-6;
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        std::fill( delta.begin( ), delta.end( ), 0. );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        error = stressTools::druckerPragerSurface( stressVector + delta, A, B, dpYield, jacobianVectorJ );

        BOOST_CHECK( ! error );

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( jacobianVectorJ - jacobianVector )/delta[ i ];

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( djacobiandstress[ j ][ i ], gradCol[ j ] ) );
        }
    }

    //Test the computation of the DP yield, jacobian, and the derivative of the jacobian
    //w.r.t. the stress for the parameter vector interface
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( jacobianVectorJ.begin( ), jacobianVectorJ.end( ), 0. );
    for ( unsigned int i=0; i<djacobiandstress.size( ); i++ ){
        std::fill( djacobiandstress[ i ].begin( ), djacobiandstress[ i ].end( ), 0. );
    }

    error = stressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, djacobiandstress );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) );

    eps = 1e-6;
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        std::fill( delta.begin( ), delta.end( ), 0. );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        error = stressTools::druckerPragerSurface( stressVector + delta, dpParam, dpYield, jacobianVectorJ );

        BOOST_CHECK( ! error );

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( jacobianVectorJ - jacobianVector )/delta[ i ];

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( djacobiandstress[ j ][ i ], gradCol[ j ] ) );
        }
    }

    //Test the computation of the DP yield, jacobian, unit direction, and the jacobian
    //of the unit direction jacobian.
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    for ( unsigned int i=0; i<unitDirectionJacobian.size( ); i++ ){
        std::fill( unitDirectionJacobian[ i ].begin( ), unitDirectionJacobian[ i ].end( ), 0. );
    }
    std::fill( unitDirectionVectorJ.begin( ), unitDirectionVectorJ.end( ), 0. );
    error = stressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) &&
                 vectorTools::fuzzyEquals( unitDirectionVector, unitDirectionVectorExpected ) );

    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        error = stressTools::druckerPragerSurface( stressVector + delta, A, B, dpYield, jacobianVector, unitDirectionVectorJ );

        BOOST_CHECK( ! error );

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( unitDirectionVectorJ - unitDirectionVector )/delta[ i ];

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( unitDirectionJacobian[ j ][ i ], gradCol[ j ] ) );
        }
    }

    //Test the computation of the DP yield, jacobian, unit direction, and the jacobian
    //of the unit direction jacobian from the parameter vector interface
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    for ( unsigned int i=0; i<unitDirectionJacobian.size( ); i++ ){
        std::fill( unitDirectionJacobian[ i ].begin( ), unitDirectionJacobian[ i ].end( ), 0. );
    }
    std::fill( unitDirectionVectorJ.begin( ), unitDirectionVectorJ.end( ), 0. );
    error = stressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( dpYield, dpYieldExpected ) &&
                 vectorTools::fuzzyEquals( jacobianVector, jacobianVectorExpected ) &&
                 vectorTools::fuzzyEquals( unitDirectionVector, unitDirectionVectorExpected ) );

    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        error = stressTools::druckerPragerSurface( stressVector + delta, dpParam, dpYield, jacobianVector, unitDirectionVectorJ );

        BOOST_CHECK( ! error );

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( unitDirectionVectorJ - unitDirectionVector )/delta[ i ];

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( unitDirectionJacobian[ j ][ i ], gradCol[ j ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testLinearViscoelasticity ){
    /*!
     * Test the implementation of linear finite-deformation
     * viscoelasticity.
     */

    floatType previousTime = 0.6;
    floatType currentTime = 23.8;

    floatVector previousStrain = {  3.03768940e-01,  4.54626930e-17, -3.71231060e-01,
                                    4.54626930e-17,  2.00000000e-01, -1.03633416e-16,
                                   -3.71231060e-01, -1.03633416e-16,  1.04623106e+00 };

    floatVector currentStrain = {  1.0163119 , -0.57654737,  0.33286978,
                                  -0.57654737,  0.45526608, -0.22243347,
                                   0.33286978, -0.22243347,  0.19842203 };

    floatVector previousStateVariables = { 1, 2, 3,
                                           2, 4, 5,
                                           3, 5, 6,
                                           1.0, 0.3, 0.2,
                                           0.3, 2.0, 0.1,
                                           0.2, 0.1, 3.0 };

    floatVector materialParameters = { 100, 1., 10, 150, 200 };

    floatType currentRateModifier = 2.7;
    floatType previousRateModifier = 3.5;

    // Compute the previous stress value
    floatVector previousStress = materialParameters[ 0 ] *   previousStrain
                               + materialParameters[ 3 ] * ( previousStrain - floatVector( previousStateVariables.begin( ),
                                                                                           previousStateVariables.begin( ) + 9 ) )
                               + materialParameters[ 4 ] * ( previousStrain - floatVector( previousStateVariables.begin( ) + 9,
                                                                                           previousStateVariables.begin( ) + 18 ) );

    floatType alpha = 0.5; //Trapezoidal rule

    floatVector stress;
    floatVector currentStateVariables;

    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                        previousTime, previousStrain,
                                        currentRateModifier, previousRateModifier,
                                        previousStateVariables,
                                        materialParameters, alpha,
                                        stress, currentStateVariables );

    //!Test for symmetry in the output stress
    BOOST_CHECK( vectorTools::fuzzyEquals( stress[ 1 ], stress[ 3 ] ) &&
                 vectorTools::fuzzyEquals( stress[ 2 ], stress[ 6 ] ) &&
                 vectorTools::fuzzyEquals( stress[ 5 ], stress[ 7 ] ) );

    //!Check that dStress is being computed correctly
    floatVector dStress;
    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                        previousTime, previousStrain,
                                        currentRateModifier, previousRateModifier,
                                        previousStateVariables,
                                        materialParameters, alpha,
                                        dStress, stress, currentStateVariables );

    BOOST_CHECK( vectorTools::fuzzyEquals( dStress, stress - previousStress ) );    

    //!Test for passing the state variables through properly
    currentTime = previousTime;

    errorOut res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                       previousTime, previousStrain,
                                                       currentRateModifier, previousRateModifier,
                                                       previousStateVariables,
                                                       materialParameters, alpha,
                                                       stress, currentStateVariables );

    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousStateVariables, currentStateVariables ) );

    //!Test to make sure the state variable evolution is occurring
    //!as expected for very large dt and alpha=0 ( fully implicit )

    currentTime = 1e10;

    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, 0.,
                                              stress, currentStateVariables );

    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( currentStrain*materialParameters[ 0 ], stress ) );

    unsigned int dim = currentStrain.size( );

    for ( unsigned int i=0; i<currentStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector subv;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        vectorTools::getValuesByIndex( currentStateVariables, indices, subv );
        BOOST_CHECK( vectorTools::fuzzyEquals( subv, currentStrain ) );
    }

    //!Test for frame invariance
    floatVector Q = { -0.44956296, -0.88488713, -0.12193405,
                      -0.37866166,  0.31242661, -0.87120891,
                       0.80901699, -0.3454915 , -0.47552826 };

    floatVector QT( Q.size( ), 0 );
    for ( unsigned int i=0; i<3; i++ ){
        for ( unsigned int j=0; j<3; j++ ){
            QT[ 3*j + i ] = Q[ 3*i + j ];
        }
    }

    //!Rotate the previous strain
    floatVector rotatedPreviousStrain( currentStrain.size( ), 0 );
    res = constitutiveTools::rotateMatrix( previousStrain, Q, rotatedPreviousStrain );

    BOOST_CHECK( ! res );

    //!Rotate the current strain
    floatVector rotatedCurrentStrain( currentStrain.size( ), 0 );
    res = constitutiveTools::rotateMatrix( currentStrain, Q, rotatedCurrentStrain );

    BOOST_CHECK( ! res );

    //!Rotate the previous state variables
    floatVector rotatedPreviousStateVariables;
    for ( unsigned int i=0; i<previousStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector subv, rotatedSubv;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        vectorTools::getValuesByIndex( previousStateVariables, indices, subv );
        res = constitutiveTools::rotateMatrix( subv, Q, rotatedSubv );
        BOOST_CHECK( ! res );
        rotatedPreviousStateVariables = vectorTools::appendVectors( {rotatedPreviousStateVariables, rotatedSubv } );
    }

    currentTime = 23.8;

    floatVector rotatedStress;
    floatVector rotatedCurrentStateVariables;

    //Calculate the rotated stress
    res = stressTools::linearViscoelasticity( currentTime,  rotatedCurrentStrain,
                                              previousTime, rotatedPreviousStrain,
                                              currentRateModifier, previousRateModifier,
                                              rotatedPreviousStateVariables,
                                              materialParameters, alpha,
                                              rotatedStress, rotatedCurrentStateVariables );

    BOOST_CHECK( ! res );

    //Re-calculate the initial values
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stress, currentStateVariables );

    BOOST_CHECK( ! res );

    floatVector stresspp;
    res = constitutiveTools::rotateMatrix( rotatedStress, QT, stresspp );
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( stress, stresspp ) );

    for ( unsigned int i=0; i<rotatedCurrentStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector rCSV, CSV, CSVpp;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        vectorTools::getValuesByIndex( rotatedCurrentStateVariables, indices, rCSV );
        vectorTools::getValuesByIndex( currentStateVariables, indices, CSV );

        res = constitutiveTools::rotateMatrix( rCSV, QT, CSVpp );
        BOOST_CHECK( ! res );
        BOOST_CHECK( vectorTools::fuzzyEquals( CSVpp, CSV ) );
    }

    //Test the implementation of the jacobian
    floatType eps = 1e-6;
    floatVector deltaStress, deltaDStress;

    floatMatrix jacobian;
    floatVector dstressdrateModifier;

    //Compute the jacobian
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stress, currentStateVariables );

    floatVector stressJ, dStressJ, currentStateVariablesJ;
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stressJ, currentStateVariablesJ, jacobian,
                                              dstressdrateModifier );

    BOOST_CHECK( vectorTools::fuzzyEquals( stressJ, stress ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( currentStateVariables, currentStateVariablesJ ) );

    for ( unsigned int i=0; i<currentStrain.size( ); i++ ){
        floatVector deltaStrain( currentStrain.size( ), 0 );
        deltaStrain[ i ] = fabs( eps*currentStrain[ i ] );

        res = stressTools::linearViscoelasticity( currentTime,  currentStrain + deltaStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  deltaStress, currentStateVariablesJ );
        BOOST_CHECK( ! res );

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<deltaStress.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( jacobian[ j ][ i ], ( deltaStress[ j ] - stressJ[ j ] )/deltaStrain[ i ] ) );
        }

    }

    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              dStressJ, stressJ, currentStateVariablesJ, jacobian,
                                              dstressdrateModifier );

    BOOST_CHECK( vectorTools::fuzzyEquals( stressJ, stress ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dStressJ, dStress ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( currentStateVariables, currentStateVariablesJ ) );

    for ( unsigned int i=0; i<currentStrain.size( ); i++ ){
        floatVector deltaStrain( currentStrain.size( ), 0 );
        deltaStrain[ i ] = fabs( eps*currentStrain[ i ] );

        res = stressTools::linearViscoelasticity( currentTime,  currentStrain + deltaStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  deltaStress, currentStateVariablesJ );
        BOOST_CHECK( ! res );

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<deltaStress.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( jacobian[ j ][ i ], ( deltaStress[ j ] - stressJ[ j ] )/deltaStrain[ i ] ) );
        }

        res = stressTools::linearViscoelasticity( currentTime,  currentStrain + deltaStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  deltaDStress, deltaStress, currentStateVariablesJ );
        BOOST_CHECK( ! res );

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<deltaStress.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( jacobian[ j ][ i ], ( deltaDStress[ j ] - dStressJ[ j ] )/deltaStrain[ i ] ) );
        }

    }

    floatVector deltaStrain( currentStrain.size( ), 0 );

    //Test the gradient w.r.t. the rate modifier
    res = stressTools::linearViscoelasticity( currentTime, currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier + fabs( eps*currentRateModifier ), previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              deltaStress, currentStateVariables );

    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( (deltaStress - stress )/fabs( eps*currentRateModifier ), dstressdrateModifier ) );

    res = stressTools::linearViscoelasticity( currentTime, currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier + fabs( eps*currentRateModifier ), previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              deltaDStress, deltaStress, currentStateVariablesJ );

    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( (deltaStress - stress )/fabs( eps*currentRateModifier ), dstressdrateModifier ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( (deltaDStress - dStress )/fabs( eps*currentRateModifier ), dstressdrateModifier ) );

    //Test to make sure odd numbers of prony series terms can be passed in
    floatVector materialParametersOdd = { materialParameters[ 1 ], 1, 10, 100, 400, 300, 200 };
    floatVector previousStateVariablesOdd( 3 * 9, 0 );

    res = stressTools::linearViscoelasticity(  currentTime, currentStrain,
                                               previousTime, previousStrain,
                                               currentRateModifier, previousRateModifier,
                                               previousStateVariablesOdd,
                                               materialParametersOdd, alpha,
                                               stress, currentStateVariables );

    BOOST_CHECK( ! res );


}

BOOST_AUTO_TEST_CASE( testVolumetricNeoHookean ){
    /*!
     * Test the computation of the mean stress ( -pressure ) using a Neo-Hookean model
     */

    //Define the bulk modulus
    floatType bulkModulus = 150.;

    //Define the deformation gradient
    floatVector deformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    //Run the test when there is no deformation
    floatType meanStress;

    errorOut res = stressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, 0. ) );

    //Test the meanStress computation as compared to the expected value
    deformationGradient = { 0.39874077,  0.11561812, -0.75485222,
                            0.14034205,  0.15851022,  1.29640525,
                            0.26235075, -0.26051883,  0.45378251 };

    floatType J = 0.25430054895115856;

    res = stressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, 0.5*bulkModulus*( J - 1/J ) ) );

    //Test the meanStress computation subject to a rotation
    deformationGradient = { -0.2350804 ,  0.16410551, -1.13402371,
                             0.1296644 , -0.22975865, -1.03460443,
                            -0.4188584 , -0.16322821,  0.31618178 };

    res = stressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals( meanStress, 0.5*bulkModulus*( J - 1/J ) ) );

    //Test the computation of the derivative of the mean stress w.r.t. J
    floatType jacobianMeanStress;
    floatType dmeanStressdJ;
    res = stressTools::volumetricNeoHookean( deformationGradient, bulkModulus, jacobianMeanStress, dmeanStressdJ );
    BOOST_CHECK( ! res );

    //Make sure the mean stress is identical
    BOOST_CHECK( vectorTools::fuzzyEquals( jacobianMeanStress, meanStress ) );

    //Make sure the jacobian is correct
    floatType ms;
    floatType eps = 1e-6;
    floatVector dmeanStressdF( deformationGradient.size( ), 0 );
    for ( unsigned int i=0; i<deformationGradient.size( ); i++ ){
        floatVector delta( deformationGradient.size( ), 0 );
        delta[ i ] = fabs( eps*deformationGradient[ i ] );
        stressTools::volumetricNeoHookean( deformationGradient + delta, bulkModulus, ms );

        dmeanStressdF[ i ] = ( ms - meanStress )/delta[ i ];
    }

    floatVector dJdF = vectorTools::computeDDetAdJ( deformationGradient, 3, 3 );

    BOOST_CHECK( vectorTools::fuzzyEquals( dmeanStressdJ*dJdF, dmeanStressdF ) );


}

BOOST_AUTO_TEST_CASE( testPeryznaModel ){
    /*!
     * Test the implementation of the Peryzna style model
     */

    floatType f = 2.;
    floatType q = 5.42;
    floatType A = 1.4;
    floatType n = 2.4;
    floatVector parameters = { n };

    floatType p;
    errorOut error = stressTools::peryznaModel( f, q, A, n, p );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( p, A*pow( (f/q ), n ) ) );

    error = stressTools::peryznaModel( f, q, A, parameters, p );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( p, A*pow( (f/q ), n ) ) );

    floatType pJ;
    floatType dpdf, dpdq, dpdA;
    error = stressTools::peryznaModel( f, q, A, n, pJ, dpdf, dpdq, dpdA );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( p, pJ ) );

    floatType eps = 1e-6;
    floatType delta = eps*fabs( f ) + eps;
    error = stressTools::peryznaModel( f + delta, q, A, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdf, 1e-5, 1e-5 ) );

    delta = eps*fabs( q ) + eps;
    error = stressTools::peryznaModel( f, q + delta, A, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdq, 1e-5, 1e-5 ) );

    delta = eps*fabs( A ) + eps;
    error = stressTools::peryznaModel( f, q, A + delta, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdA, 1e-5, 1e-5 ) );

    f = -1;
    error = stressTools::peryznaModel( f, q, A, n, p );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( p, 0. ) );

    error = stressTools::peryznaModel( f, q, A, n, pJ, dpdf, dpdq, dpdA );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( p, pJ ) );

    delta = eps*fabs( f ) + eps;
    error = stressTools::peryznaModel( f + delta, q, A, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdf, 1e-5, 1e-5 ) );

    delta = eps*fabs( q ) + eps;
    error = stressTools::peryznaModel( f, q + delta, A, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdq, 1e-5, 1e-5 ) );

    delta = eps*fabs( A ) + eps;
    error = stressTools::peryznaModel( f, q, A + delta, n, pJ );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( (pJ - p )/delta, dpdA, 1e-5, 1e-5 ) );

    floatType pJv, dpdfv, dpdqv, dpdAv;
    error = stressTools::peryznaModel( f, q, A, parameters, pJv, dpdfv, dpdqv, dpdAv );

    BOOST_CHECK( vectorTools::fuzzyEquals( pJv, p ) &&
                 vectorTools::fuzzyEquals( dpdfv, dpdf ) &&
                 vectorTools::fuzzyEquals( dpdqv, dpdq ) &&
                 vectorTools::fuzzyEquals( dpdAv, dpdA ) );

}

BOOST_AUTO_TEST_CASE( testLinearHardening ){
    /*!
     * Test the linear hardening function.
     */

    floatVector stateVariables = { 1, 2, 3, 4, 5 };
    floatVector linearModuli = { 6, 7, 8, 9, 10 };
    floatType shiftFactor = 3.7;
    floatType value;

    errorOut error = stressTools::linearHardening( stateVariables, linearModuli, shiftFactor, value );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( value, 133.7 ) );

    floatType valueJ;
    floatVector valueJacobian;

    error = stressTools::linearHardening( stateVariables, linearModuli, shiftFactor, valueJ, valueJacobian );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( value, valueJ ) );

    floatType eps = 1e-6;
    for ( unsigned int i=0; i<stateVariables.size( ); i++ ){
        floatVector delta( stateVariables.size( ), 0 );
        delta[ i ] = eps*fabs( stateVariables[ i ] ) + eps;

        error = stressTools::linearHardening( stateVariables + delta, linearModuli, shiftFactor, valueJ );

        BOOST_CHECK( ! error );

        BOOST_CHECK( vectorTools::fuzzyEquals( (valueJ - value )/delta[ i ], valueJacobian[ i ] ) );
    }

}

BOOST_AUTO_TEST_CASE( testComputeJaumannStiffnessTensor ){
    /*!
     * Test the computation of the Jaumann stiffness tensor
     */

    floatVector cauchyStress = { 0.39293837, -0.42772133, -0.54629709,
                                 0.10262954,  0.43893794, -0.15378708,
                                 0.9615284 ,  0.36965948, -0.0381362 };

    floatVector currentDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                              -0.06142776,  0.5596779 , -0.10195574,
                                               0.23799541, -0.31750827,  0.67545176 };

    floatMatrix dCauchydF =
    {
        { 0.3155137384183837, 0.3182758709686606, 1.3440095855132106, 3.4943179407778957, 2.244553248606352, 1.1102351067758287, 2.2244338257022154, -1.770410861468218, -1.382113443776859 },
        { -2.717367691210444, -2.0628595361117066, 1.3097612385448776, -4.078950600549248, -0.6629882732047176, -0.6913723667035621, -0.06314902349693785, -0.7416970970417203, -1.8773877702753472 },
        { -0.7364869303719179, 3.933891631171348, 4.441600182038796, 0.018366758843365583, 1.2395295179211119, -3.843816049207043, -1.827145181796791, -0.8517378804636822, 3.6630915788336593 },
        { -2.495446346034933, -0.1696573573729565, 4.85559785610705, 0.19485119259809336, 1.1289452576296766, -3.7937133400967626, 3.2634080050683325, 1.0306012841092738, 0.45068006466464916 },
        { -1.5723616622569159, -1.958792109728159, -0.8297778897529839, 1.8130076579279664, 3.754568417951749, 0.10422337478011134, 1.6931378296227229, 0.8593655256221289, 1.2490350209559986 },
        { 1.7468905098782483, 3.423424376202573, -4.168050116675612, 2.6368284144333822, -2.5633362546312597, -3.0577703942122914, 0.7245695749147307, -4.0428748338761284, 3.8532682627513957 },
        { 1.2724897205126873, 2.234163581899548, -4.8387079330498315, 0.944318794450425, 0.567851923942887, -3.4104035585527726, -3.469294848752269, 1.955295287709109, -1.8123357361812364 },
        { 1.9197029553181966, 0.5438324971777209, -1.110494258768554, 4.251324896139861, 3.416699969127163, -1.4260243331682376, -4.564085362009594, -1.9523192658890254, -1.0181431808201902 },
        { 2.0495883045136223, 4.953584820340175, -1.4408513428254044, 2.625478137854338, 0.9317691656222116, 1.9170179870017712, -3.488725476519198, -1.0112370727384346, -2.591441022763755 }
    };

    floatMatrix CAnswer =
    {
        { 0.9323460532562229, 1.5206835766443023, 1.413686140800323, 1.5206835766443023, 1.3213222576565693, -0.05885487383212795, 1.413686140800323, -0.05885487383212795, 0.5509125382847507 },
        { -2.2284304975160723, -2.4072217797743627, 0.2614487023857131, -2.4072217797743627, -0.47773104837760527, -0.7235389116516223, 0.2614487023857131, -0.7235389116516223, -1.4753404121609062 },
        { -0.8029034532123083, 0.36784058323148006, 0.45916803313893, 0.36784058323148006, 0.5382110794054509, -1.8617183183800223, 0.45916803313893, -1.8617183183800223, 1.7635262086408587 },
        { -0.9848225913718358, -0.6544410518218504, 2.796317798370991, -0.6544410518218504, 1.1092968426129497, -1.2720772977623038, 2.796317798370991, -1.2720772977623038, 0.8564938612296599 },
        { -0.8466721625112175, 0.06869199010472715, 0.6745222632776061, 0.06869199010472715, 2.418291731814228, -0.2203026144049039, 0.6745222632776061, -0.2203026144049039, 1.4127042017386933 },
        { -0.08691431534171057, 2.1438154641991627, -0.6617637870366262, 2.1438154641991627, -1.4386469161393038, -1.6620110384171627, -0.6617637870366262, -1.6620110384171627, 3.9050001576917728 },
        { 0.6380681668217065, 0.8189097972242327, -3.7457997213743823, 0.8189097972242327, 1.5690454177333444, -0.38344631567380594, -3.7457997213743823, -0.38344631567380594, -1.7091156173490465 },
        { 1.74261683722654, 1.6149532727880505, -2.2323077259418067, 1.6149532727880505, 2.166152956534386, -0.8723741082663425, -2.2323077259418067, -0.8723741082663425, -0.7844009581372655 },
        { 0.6834798573624237, 2.714309478004904, -2.8027874159759896, 2.714309478004904, 0.12662618464897926, 0.7682051607516402, -2.8027874159759896, 0.7682051607516402, -2.2977540890151382 }
    };

    floatMatrix CResult;

    errorOut error = stressTools::computeJaumannStiffnessTensor( cauchyStress, currentDeformationGradient, dCauchydF, CResult );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( CResult, CAnswer ) );

}
