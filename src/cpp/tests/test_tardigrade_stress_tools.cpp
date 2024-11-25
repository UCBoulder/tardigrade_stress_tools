/**
  * \file test_tardigrade_stress_tools.cpp
  *
  * Tests for tardigrade_stress_tools
  */

#include<tardigrade_stress_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_stress_tools
#include <boost/test/included/unit_test.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeConstitutiveTools::floatType floatType;
typedef tardigradeConstitutiveTools::floatVector floatVector;
typedef tardigradeConstitutiveTools::floatMatrix floatMatrix;

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

BOOST_AUTO_TEST_CASE( testCalculateMeanStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
    tardigradeStressTools::calculateMeanStress( stressVector, meanStress );
    

    BOOST_TEST( meanStress == meanStressExpected );

    //Test for correct mean stress calculation from stressVector without pointer output
    meanStress = 0.;
    meanStress = tardigradeStressTools::calculateMeanStress( stressVector );
    BOOST_TEST( meanStress == meanStressExpected );

    //Test for correct mean stress calculation from stressMatrix with pointer output
    meanStress = 0.;
    tardigradeStressTools::calculateMeanStress( stressMatrix, meanStress );
    BOOST_TEST( meanStress == meanStressExpected );

    //Test for correct mean stress calculation from stressMatrix without pointer output
    meanStress = 0.;
    meanStress = tardigradeStressTools::calculateMeanStress( stressMatrix );
    BOOST_TEST( meanStress == meanStressExpected );

    //Test for correct mean stress and jacobian calculation
    meanStress = 0.;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    tardigradeStressTools::calculateMeanStress( stressVector, meanStress, jacobianVector );
    BOOST_TEST( meanStress == meanStressExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( testCalculateDeviatoricStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
    

    //Test computation of deviatoric tensor in row major format
    std::fill( deviatoricVector.begin( ), deviatoricVector.end( ), 0. );
    tardigradeStressTools::calculateDeviatoricStress( stressVector, deviatoricVector );

    

    BOOST_TEST( expectedVector == deviatoricVector, CHECK_PER_ELEMENT );

    std::fill( deviatoricVector.begin( ), deviatoricVector.end( ), 0. );

    deviatoricVector = tardigradeStressTools::calculateDeviatoricStress( stressVector );

    BOOST_TEST( expectedVector == deviatoricVector, CHECK_PER_ELEMENT );

    //Test the computation of the jacobian
    tardigradeStressTools::calculateDeviatoricStress( stressVector, deviatoricVectorJ, jacobian );

    

    BOOST_TEST( expectedVector == deviatoricVectorJ, CHECK_PER_ELEMENT );

    //Test the jacobian
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps * fabs( stressVector[ i ] ) + eps;

        floatVector Rp( 9, 0 ), Rm( 9, 0 );

        tardigradeStressTools::calculateDeviatoricStress( stressVector + delta, Rp );

        

        tardigradeStressTools::calculateDeviatoricStress( stressVector - delta, Rm );

        

        floatVector grad = ( Rp - Rm ) / ( 2 * delta[ i ] );

        for ( unsigned int j=0; j<grad.size( ); j++ ){
            BOOST_TEST( grad[ j ] == jacobian[ j ][ i ] );
        }

    }

    deviatoricVector = tardigradeStressTools::calculateDeviatoricStress( stressVector, jacobian2 );

    BOOST_TEST( deviatoricVector == expectedVector, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian ) == tardigradeVectorTools::appendVectors( jacobian2 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( testCalculateVonMisesStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
    tardigradeStressTools::calculateVonMisesStress( stressVector, vonMises );
    BOOST_TEST( vonMises == vonMisesExpected );

    vonMises = 0.;
    vonMises = tardigradeStressTools::calculateVonMisesStress( stressVector );
    BOOST_TEST( vonMises == vonMisesExpected );

    //Test computation of vonMises stress and jacobian
    vonMises = 0.;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    tardigradeStressTools::calculateVonMisesStress( stressVector, vonMises, jacobianVector );
    BOOST_TEST( vonMises == vonMisesExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( testDruckerPragerSurface, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
    tardigradeStressTools::druckerPragerSurface( vonMises, meanStress, A, B, dpYield );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    tardigradeStressTools::druckerPragerSurface( vonMises, meanStress, dpParam, dpYield );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    dpYield = tardigradeStressTools::druckerPragerSurface( vonMises, meanStress, A, B );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    dpYield = tardigradeStressTools::druckerPragerSurface( vonMises, meanStress, dpParam );
    BOOST_TEST( dpYield == dpYieldExpected );

    //Test computation of DP yield criterion from row major stress tensor
    dpYield = 0;
    tardigradeStressTools::druckerPragerSurface( stressVector, A, B, dpYield );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    tardigradeStressTools::druckerPragerSurface( stressVector, dpParam, dpYield );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    dpYield = tardigradeStressTools::druckerPragerSurface( stressVector, A, B );
    BOOST_TEST( dpYield == dpYieldExpected );

    dpYield = 0;
    dpYield = tardigradeStressTools::druckerPragerSurface( stressVector, dpParam );
    BOOST_TEST( dpYield == dpYieldExpected );

    //Test computation of DP yield and jacobian from row major stress tensor
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    tardigradeStressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector );
    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );

    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    tardigradeStressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector );
    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );

    //Test computation of DP yield, jacobian, and unit normal from row major stress tensor
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( unitDirectionVector.begin( ), unitDirectionVector.end( ), 0. );
    tardigradeStressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, unitDirectionVector );
    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );
    BOOST_TEST( unitDirectionVector == unitDirectionVectorExpected, CHECK_PER_ELEMENT );

    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( unitDirectionVector.begin( ), unitDirectionVector.end( ), 0. );
    tardigradeStressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector );
    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );
    BOOST_TEST( unitDirectionVector == unitDirectionVectorExpected, CHECK_PER_ELEMENT );

    //Test the computation of the DP yield, jacobian, and the derivative of the jacobian
    //w.r.t. the stress
    dpYield = 0;
    std::fill( jacobianVector.begin( ), jacobianVector.end( ), 0. );
    std::fill( jacobianVectorJ.begin( ), jacobianVectorJ.end( ), 0. );
    for ( unsigned int i=0; i<djacobiandstress.size( ); i++ ){
        std::fill( djacobiandstress[ i ].begin( ), djacobiandstress[ i ].end( ), 0. );
    }

    tardigradeStressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, djacobiandstress );

    

    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );

    eps = 1e-6;
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        std::fill( delta.begin( ), delta.end( ), 0. );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        floatVector jvp(9,0), jvm(9,0);

        tardigradeStressTools::druckerPragerSurface( stressVector + delta, A, B, dpYield, jvp );

        

        tardigradeStressTools::druckerPragerSurface( stressVector - delta, A, B, dpYield, jvm );

        

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( jvp - jvm )/(2*delta[ i ]);

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_TEST( djacobiandstress[ j ][ i ] == gradCol[ j ] );
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

    tardigradeStressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, djacobiandstress );

    

    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected );

    eps = 1e-6;
    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        std::fill( delta.begin( ), delta.end( ), 0. );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        floatVector jvp(9,0), jvm(9,0);

        tardigradeStressTools::druckerPragerSurface( stressVector + delta, dpParam, dpYield, jvp );

        

        tardigradeStressTools::druckerPragerSurface( stressVector - delta, dpParam, dpYield, jvm );

        

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( jvp - jvm )/(2*delta[ i ]);

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_TEST( djacobiandstress[ j ][ i ] == gradCol[ j ] );
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
    tardigradeStressTools::druckerPragerSurface( stressVector, A, B, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian );

    

    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected );
    BOOST_TEST( unitDirectionVector == unitDirectionVectorExpected, CHECK_PER_ELEMENT );

    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        floatVector uvp(9,0), uvm(9,0);

        tardigradeStressTools::druckerPragerSurface( stressVector + delta, A, B, dpYield, jacobianVector, uvp );

        

        tardigradeStressTools::druckerPragerSurface( stressVector - delta, A, B, dpYield, jacobianVector, uvm );

        

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( uvp - uvm ) / ( 2 * delta[ i ] );

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_TEST( unitDirectionJacobian[ j ][ i ] == gradCol[ j ] );
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
    tardigradeStressTools::druckerPragerSurface( stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian );

    

    BOOST_TEST( dpYield == dpYieldExpected );
    BOOST_TEST( jacobianVector == jacobianVectorExpected, CHECK_PER_ELEMENT );
    BOOST_TEST( unitDirectionVector == unitDirectionVectorExpected, CHECK_PER_ELEMENT );

    for ( unsigned int i=0; i<stressVector.size( ); i++ ){
        floatVector delta( stressVector.size( ), 0 );
        delta[ i ] = eps*fabs( stressVector[ i ] ) + eps;

        floatVector uvp(9,0), uvm(9,0);

        tardigradeStressTools::druckerPragerSurface( stressVector + delta, dpParam, dpYield, jacobianVector, uvp );

        

        tardigradeStressTools::druckerPragerSurface( stressVector - delta, dpParam, dpYield, jacobianVector, uvm );

        

        std::fill( gradCol.begin( ), gradCol.end( ), 0. );
        gradCol = ( uvp - uvm ) / ( 2 * delta[ i ] );

        for ( unsigned int j=0; j<gradCol.size( ); j++ ){
            BOOST_TEST( unitDirectionJacobian[ j ][ i ] == gradCol[ j ] );
        }
    }

}

BOOST_AUTO_TEST_CASE( testLinearViscoelasticity, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                        previousTime, previousStrain,
                                        currentRateModifier, previousRateModifier,
                                        previousStateVariables,
                                        materialParameters, alpha,
                                        stress, currentStateVariables );

    //!Test for symmetry in the output stress
    BOOST_TEST( stress[ 1 ] == stress[ 3 ] );
    BOOST_TEST( stress[ 2 ] == stress[ 6 ] );
    BOOST_TEST( stress[ 5 ] == stress[ 7 ] );

    //!Check that dStress is being computed correctly
    floatVector dStress;
    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                        previousTime, previousStrain,
                                        currentRateModifier, previousRateModifier,
                                        previousStateVariables,
                                        materialParameters, alpha,
                                        dStress, stress, currentStateVariables );

    BOOST_TEST( dStress == (stress - previousStress ), CHECK_PER_ELEMENT );

    //!Test for passing the state variables through properly
    currentTime = previousTime;

    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                        previousTime, previousStrain,
                                        currentRateModifier, previousRateModifier,
                                        previousStateVariables,
                                        materialParameters, alpha,
                                        stress, currentStateVariables );

    BOOST_TEST( previousStateVariables == currentStateVariables, CHECK_PER_ELEMENT );

    //!Test to make sure the state variable evolution is occurring
    //!as expected for very large dt and alpha=0 ( fully implicit )

    currentTime = 1e10;

    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, 0.,
                                              stress, currentStateVariables );

    

    BOOST_TEST( ( currentStrain * materialParameters[ 0 ] ) == stress, CHECK_PER_ELEMENT );

    unsigned int dim = currentStrain.size( );

    for ( unsigned int i=0; i<currentStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector subv;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        tardigradeVectorTools::getValuesByIndex( currentStateVariables, indices, subv );
        BOOST_TEST( subv == currentStrain, CHECK_PER_ELEMENT );
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
    tardigradeConstitutiveTools::rotateMatrix( previousStrain, Q, rotatedPreviousStrain );

    

    //!Rotate the current strain
    floatVector rotatedCurrentStrain( currentStrain.size( ), 0 );
    tardigradeConstitutiveTools::rotateMatrix( currentStrain, Q, rotatedCurrentStrain );

    

    //!Rotate the previous state variables
    floatVector rotatedPreviousStateVariables;
    for ( unsigned int i=0; i<previousStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector subv, rotatedSubv;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        tardigradeVectorTools::getValuesByIndex( previousStateVariables, indices, subv );
        tardigradeConstitutiveTools::rotateMatrix( subv, Q, rotatedSubv );
        
        rotatedPreviousStateVariables = tardigradeVectorTools::appendVectors( {rotatedPreviousStateVariables, rotatedSubv } );
    }

    currentTime = 23.8;

    floatVector rotatedStress;
    floatVector rotatedCurrentStateVariables;

    //Calculate the rotated stress
    tardigradeStressTools::linearViscoelasticity( currentTime,  rotatedCurrentStrain,
                                              previousTime, rotatedPreviousStrain,
                                              currentRateModifier, previousRateModifier,
                                              rotatedPreviousStateVariables,
                                              materialParameters, alpha,
                                              rotatedStress, rotatedCurrentStateVariables );

    

    //Re-calculate the initial values
    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stress, currentStateVariables );

    

    floatVector stresspp;
    tardigradeConstitutiveTools::rotateMatrix( rotatedStress, QT, stresspp );
    

    BOOST_TEST( stress == stresspp, CHECK_PER_ELEMENT );

    for ( unsigned int i=0; i<rotatedCurrentStateVariables.size( )/dim; i++ ){
        std::vector< unsigned int > indices( dim, 0 );
        floatVector rCSV, CSV, CSVpp;

        for ( unsigned int j=0; j<dim; j++ ){
            indices[ j ] = i*dim + j;
        }

        tardigradeVectorTools::getValuesByIndex( rotatedCurrentStateVariables, indices, rCSV );
        tardigradeVectorTools::getValuesByIndex( currentStateVariables, indices, CSV );

        tardigradeConstitutiveTools::rotateMatrix( rCSV, QT, CSVpp );
        
        BOOST_TEST( CSVpp == CSV, CHECK_PER_ELEMENT );
    }

    //Test the implementation of the jacobian
    floatType eps = 1e-6;
    floatVector deltaStress, deltaDStress;

    floatMatrix jacobian;
    floatVector dstressdrateModifier;

    //Compute the jacobian
    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stress, currentStateVariables );

    

    floatVector stressJ, dStressJ, currentStateVariablesJ;
    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              stressJ, currentStateVariablesJ, jacobian,
                                              dstressdrateModifier );

    
    BOOST_TEST( stressJ == stress, CHECK_PER_ELEMENT );

    BOOST_TEST( currentStateVariables == currentStateVariablesJ, CHECK_PER_ELEMENT );

    floatVector stress_2, dStress_2, currentStateVariables_2, dstressdrateModifier_2;
    floatMatrix jacobian_2;

    tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables,
                                              materialParameters, alpha,
                                              dStress_2, stress_2, currentStateVariables_2, jacobian_2,
                                              dstressdrateModifier_2 );

    

    BOOST_TEST( stress == stress_2, CHECK_PER_ELEMENT );
    BOOST_TEST( dStress == dStress_2, CHECK_PER_ELEMENT );
    BOOST_TEST( currentStateVariables == currentStateVariables_2, CHECK_PER_ELEMENT );
    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian ) == tardigradeVectorTools::appendVectors( jacobian_2 ), CHECK_PER_ELEMENT );
    BOOST_TEST( dstressdrateModifier == dstressdrateModifier_2, CHECK_PER_ELEMENT );

    floatVector stress_3, dStress_3, currentStateVariables_3, dstressdrateModifier_3, dstressdPreviousRateModifier,
                dStateVariablesdRateModifier, dStateVariablesdPreviousRateModifier;

    floatMatrix jacobian_3, dstressdPreviousStrain, dstressdPreviousStateVariables,
                dStateVariablesdStrain, dStateVariablesdPreviousStrain, dStateVariablesdPreviousStateVariables;

    tardigradeStressTools::linearViscoelasticity( currentTime, currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariables, materialParameters,
                                              alpha,
                                              dStress_3, stress_3, currentStateVariables_3,
                                              jacobian_3, dstressdrateModifier_3,
                                              dstressdPreviousStrain, dstressdPreviousRateModifier,
                                              dstressdPreviousStateVariables,
                                              dStateVariablesdStrain, dStateVariablesdRateModifier,
                                              dStateVariablesdPreviousStrain, dStateVariablesdPreviousRateModifier,
                                              dStateVariablesdPreviousStateVariables );

    
    BOOST_TEST( stress_3 == stress, CHECK_PER_ELEMENT );

    BOOST_TEST( currentStateVariables_3 == currentStateVariables, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian_3 ) == tardigradeVectorTools::appendVectors( jacobian ), CHECK_PER_ELEMENT );

    BOOST_TEST( dstressdrateModifier_3 == dstressdrateModifier, CHECK_PER_ELEMENT );

    floatMatrix jacobian_answer( stress.size( ), floatVector( currentStrain.size( ), 0 ) );
    floatVector dstressdrateModifier_answer( stress.size( ), 0 );
    floatMatrix dstressdPreviousStrain_answer( stress.size( ), floatVector( currentStrain.size( ), 0 ) );
    floatVector dstressdPreviousRateModifier_answer( stress.size( ), 0 );
    floatMatrix dstressdPreviousStateVariables_answer( stress.size( ), floatVector( previousStateVariables.size( ), 0 ) );
    floatMatrix dStateVariablesdStrain_answer( currentStateVariables.size( ), floatVector( currentStrain.size( ), 0 ) );
    floatVector dStateVariablesdRateModifier_answer( currentStateVariables.size( ), 0 );
    floatMatrix dStateVariablesdPreviousStrain_answer( currentStateVariables.size( ), floatVector( previousStrain.size( ), 0 ) );
    floatVector dStateVariablesdPreviousRateModifier_answer( currentStateVariables.size( ), 0 );
    floatMatrix dStateVariablesdPreviousStateVariables_answer( currentStateVariables.size( ), floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i=0; i<currentStrain.size( ); i++ ){

        floatVector delta( currentStrain.size( ), 0 );

        delta[ i ] = eps * fabs( currentStrain[ i ] ) + eps;

        floatVector s_p, s_m, xi_p, xi_m;

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain + delta,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_p, xi_p );

        

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain - delta,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_m, xi_m );

        

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<stress.size( ); j++ ){

            jacobian_answer[ j ][ i ] = ( s_p[ j ] - s_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j=0; j<currentStateVariables.size( ); j++ ){

            dStateVariablesdStrain_answer[ j ][ i ] = ( xi_p[ j ] - xi_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian ) == tardigradeVectorTools::appendVectors( jacobian_answer ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariablesdStrain ) == tardigradeVectorTools::appendVectors( dStateVariablesdStrain_answer ), CHECK_PER_ELEMENT );

    //Test the gradient w.r.t. the rate modifier
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( currentRateModifier ) + eps;

        floatVector s_p, s_m;
        floatVector xi_p, xi_m;
    
        tardigradeStressTools::linearViscoelasticity( currentTime, currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier + delta, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_p, xi_p );

        

        tardigradeStressTools::linearViscoelasticity( currentTime, currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier - delta, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_m, xi_m );

        

        dstressdrateModifier_answer = ( s_p - s_m ) / ( 2 * delta );

        dStateVariablesdRateModifier_answer = ( xi_p - xi_m ) / ( 2 * delta );

    }

    BOOST_TEST( dstressdrateModifier == dstressdrateModifier_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( dStateVariablesdRateModifier == dStateVariablesdRateModifier_answer, CHECK_PER_ELEMENT );

    for ( unsigned int i=0; i<previousStrain.size( ); i++ ){

        floatVector delta( previousStrain.size( ), 0 );

        delta[ i ] = eps * fabs( previousStrain[ i ] ) + eps;

        floatVector s_p, s_m, xi_p, xi_m;

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                  previousTime, previousStrain + delta,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_p, xi_p );

        

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                  previousTime, previousStrain - delta,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_m, xi_m );

        

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<stress.size( ); j++ ){

            dstressdPreviousStrain_answer[ j ][ i ] = ( s_p[ j ] - s_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j=0; j<currentStateVariables.size( ); j++ ){

            dStateVariablesdPreviousStrain_answer[ j ][ i ] = ( xi_p[ j ] - xi_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dstressdPreviousStrain ) == tardigradeVectorTools::appendVectors( dstressdPreviousStrain_answer ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariablesdPreviousStrain ) == tardigradeVectorTools::appendVectors( dStateVariablesdPreviousStrain_answer ), CHECK_PER_ELEMENT );

    //Test the gradient w.r.t. the previous rate modifier
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousRateModifier ) + eps;

        floatVector s_p, s_m;
        floatVector xi_p, xi_m;
    
        tardigradeStressTools::linearViscoelasticity( currentTime, currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier + delta,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_p, xi_p );

        

        tardigradeStressTools::linearViscoelasticity( currentTime, currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier - delta,
                                                  previousStateVariables,
                                                  materialParameters, alpha,
                                                  s_m, xi_m );

        

        dstressdPreviousRateModifier_answer = ( s_p - s_m ) / ( 2 * delta );

        dStateVariablesdPreviousRateModifier_answer = ( xi_p - xi_m ) / ( 2 * delta );

    }

    BOOST_TEST( dstressdPreviousRateModifier == dstressdPreviousRateModifier_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( dStateVariablesdPreviousRateModifier == dStateVariablesdPreviousRateModifier_answer, CHECK_PER_ELEMENT );

    for ( unsigned int i=0; i<previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * fabs( previousStateVariables[ i ] ) + eps;

        floatVector s_p, s_m, xi_p, xi_m;

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables + delta,
                                                  materialParameters, alpha,
                                                  s_p, xi_p );

        

        tardigradeStressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                  previousTime, previousStrain,
                                                  currentRateModifier, previousRateModifier,
                                                  previousStateVariables - delta,
                                                  materialParameters, alpha,
                                                  s_m, xi_m );

        

        //Compare the values in the column to the jacobian's values
        for ( unsigned int j=0; j<stress.size( ); j++ ){

            dstressdPreviousStateVariables_answer[ j ][ i ] = ( s_p[ j ] - s_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j=0; j<currentStateVariables.size( ); j++ ){

            dStateVariablesdPreviousStateVariables_answer[ j ][ i ] = ( xi_p[ j ] - xi_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dstressdPreviousStateVariables ) == tardigradeVectorTools::appendVectors( dstressdPreviousStateVariables_answer ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariablesdPreviousStateVariables ) == tardigradeVectorTools::appendVectors( dStateVariablesdPreviousStateVariables_answer ), CHECK_PER_ELEMENT );

    //Test to make sure odd numbers of prony series terms can be passed in
    floatVector materialParametersOdd = { materialParameters[ 1 ], 1, 10, 100, 400, 300, 200 };
    floatVector previousStateVariablesOdd( 3 * 9, 0 );

    tardigradeStressTools::linearViscoelasticity(  currentTime, currentStrain,
                                               previousTime, previousStrain,
                                               currentRateModifier, previousRateModifier,
                                               previousStateVariablesOdd,
                                               materialParametersOdd, alpha,
                                               stress, currentStateVariables );

    


}

BOOST_AUTO_TEST_CASE( testVolumetricNeoHookean, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the computation of the mean stress ( -pressure ) using a Neo-Hookean model
     */

    //Define the bulk modulus
    floatType bulkModulus = 150.;

    //Define the deformation gradient
    floatVector deformationGradient = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    //Run the test when there is no deformation
    floatType meanStress;

    tardigradeStressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    

    BOOST_TEST( meanStress == 0. );

    //Test the meanStress computation as compared to the expected value
    deformationGradient = { 0.39874077,  0.11561812, -0.75485222,
                            0.14034205,  0.15851022,  1.29640525,
                            0.26235075, -0.26051883,  0.45378251 };

    floatType J = 0.25430054895115856;

    tardigradeStressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    

    BOOST_TEST( meanStress == 0.5 * bulkModulus * ( J - 1/J ) );

    //Test the meanStress computation subject to a rotation
    deformationGradient = { -0.2350804 ,  0.16410551, -1.13402371,
                             0.1296644 , -0.22975865, -1.03460443,
                            -0.4188584 , -0.16322821,  0.31618178 };

    tardigradeStressTools::volumetricNeoHookean( deformationGradient, bulkModulus, meanStress );
    

    BOOST_TEST( meanStress == 0.5 * bulkModulus * ( J - 1/J ) );

    //Test the computation of the derivative of the mean stress w.r.t. J
    floatType jacobianMeanStress;
    floatType dmeanStressdJ;
    tardigradeStressTools::volumetricNeoHookean( deformationGradient, bulkModulus, jacobianMeanStress, dmeanStressdJ );
    

    //Make sure the mean stress is identical
    BOOST_TEST( jacobianMeanStress == meanStress );

    //Make sure the jacobian is correct
    floatType msp, msm;
    floatType eps = 1e-6;
    floatVector dmeanStressdF( deformationGradient.size( ), 0 );
    for ( unsigned int i=0; i<deformationGradient.size( ); i++ ){
        floatVector delta( deformationGradient.size( ), 0 );
        delta[ i ] = fabs( eps*deformationGradient[ i ] );

        tardigradeStressTools::volumetricNeoHookean( deformationGradient + delta, bulkModulus, msp );

        tardigradeStressTools::volumetricNeoHookean( deformationGradient - delta, bulkModulus, msm );

        dmeanStressdF[ i ] = ( msp - msm ) / ( 2 * delta[ i ] );
    }

    floatVector dJdF = tardigradeVectorTools::computeDDetADA( deformationGradient, 3, 3 );

    BOOST_TEST( ( dmeanStressdJ * dJdF ) == dmeanStressdF, CHECK_PER_ELEMENT );


}

BOOST_AUTO_TEST_CASE( testPeryznaModel, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the implementation of the Peryzna style model
     */

    floatType f = 2.;
    floatType q = 5.42;
    floatType A = 1.4;
    floatType n = 2.4;
    floatVector parameters = { n };

    floatType p;
    tardigradeStressTools::peryznaModel( f, q, A, n, p );

    

    BOOST_TEST( p == A*pow( (f/q ), n ) );

    tardigradeStressTools::peryznaModel( f, q, A, parameters, p );

    

    BOOST_TEST( p == A*pow( (f/q ), n ) );

    floatType pJ;
    floatType dpdf, dpdq, dpdA;
    tardigradeStressTools::peryznaModel( f, q, A, n, pJ, dpdf, dpdq, dpdA );

    

    BOOST_TEST( p == pJ );

    floatType eps = 1e-6;
    floatType delta = eps*fabs( f ) + eps;
    floatType pp, pm;
    tardigradeStressTools::peryznaModel( f + delta, q, A, n, pp );

    

    tardigradeStressTools::peryznaModel( f - delta, q, A, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdf );

    delta = eps*fabs( q ) + eps;
    tardigradeStressTools::peryznaModel( f, q + delta, A, n, pp );

    

    tardigradeStressTools::peryznaModel( f, q - delta, A, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdq );

    delta = eps*fabs( A ) + eps;
    tardigradeStressTools::peryznaModel( f, q, A + delta, n, pp );

    

    tardigradeStressTools::peryznaModel( f, q, A - delta, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdA );

    f = -1;
    tardigradeStressTools::peryznaModel( f, q, A, n, p );

    

    BOOST_TEST( p == 0. );

    tardigradeStressTools::peryznaModel( f, q, A, n, pJ, dpdf, dpdq, dpdA );

    

    BOOST_TEST( p == pJ );

    delta = eps*fabs( f ) + eps;
    tardigradeStressTools::peryznaModel( f + delta, q, A, n, pp );

    

    tardigradeStressTools::peryznaModel( f - delta, q, A, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdf );

    delta = eps*fabs( q ) + eps;
    tardigradeStressTools::peryznaModel( f, q + delta, A, n, pp );

    

    tardigradeStressTools::peryznaModel( f, q - delta, A, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdq );

    delta = eps*fabs( A ) + eps;
    tardigradeStressTools::peryznaModel( f, q, A + delta, n, pp );

    

    tardigradeStressTools::peryznaModel( f, q, A - delta, n, pm );

    

    BOOST_TEST( ( pp - pm ) / ( 2 * delta ) == dpdA );

    floatType pJv, dpdfv, dpdqv, dpdAv;
    tardigradeStressTools::peryznaModel( f, q, A, parameters, pJv, dpdfv, dpdqv, dpdAv );

    BOOST_TEST( pJv == p );
    BOOST_TEST( dpdfv == dpdf );
    BOOST_TEST( dpdqv == dpdq );
    BOOST_TEST( dpdAv == dpdA );

}

BOOST_AUTO_TEST_CASE( testLinearHardening, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the linear hardening function.
     */

    floatVector stateVariables = { 1, 2, 3, 4, 5 };
    floatVector linearModuli = { 6, 7, 8, 9, 10 };
    floatType shiftFactor = 3.7;
    floatType value;

    tardigradeStressTools::linearHardening( stateVariables, linearModuli, shiftFactor, value );

    

    BOOST_TEST( value = 133.7 );

    floatType valueJ;
    floatVector valueJacobian;

    tardigradeStressTools::linearHardening( stateVariables, linearModuli, shiftFactor, valueJ, valueJacobian );

    

    BOOST_CHECK( value == valueJ );

    floatType eps = 1e-6;
    for ( unsigned int i=0; i<stateVariables.size( ); i++ ){
        floatVector delta( stateVariables.size( ), 0 );
        delta[ i ] = eps*fabs( stateVariables[ i ] ) + eps;

        floatType vp, vm;

        tardigradeStressTools::linearHardening( stateVariables + delta, linearModuli, shiftFactor, vp );

        

        tardigradeStressTools::linearHardening( stateVariables - delta, linearModuli, shiftFactor, vm );

        

        BOOST_TEST( ( vp - vm ) / ( 2 * delta[ i ] ) == valueJacobian[ i ] );

    }

}

BOOST_AUTO_TEST_CASE( testComputeJaumannStiffnessTensor, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    tardigradeStressTools::computeJaumannStiffnessTensor( cauchyStress, currentDeformationGradient, dCauchydF, CResult );

    

    BOOST_TEST( tardigradeVectorTools::appendVectors( CResult ) == tardigradeVectorTools::appendVectors( CAnswer ), CHECK_PER_ELEMENT );

}
