/**
  * \file test_linear_elasticity.cpp
  *
  * Tests for the linear elasticity support module
  */

#include<linear_elasticity.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_linear_elasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef stressTools::linearElasticity::floatType floatType; //!< Redefinition for the float type
typedef stressTools::linearElasticity::floatVector floatVector; //1< Redefinition for the float vector
typedef stressTools::linearElasticity::floatMatrix floatMatrix; //1< Redefinition for the float matrix

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( formReferenceStiffnessTensor ){

    // Unique component values to check value and location
    floatType C1111 =  0.;
    floatType C1112 =  1.;
    floatType C1113 =  2.;
    floatType C1122 =  3.;
    floatType C1123 =  4.;
    floatType C1133 =  5.;
    floatType C1212 =  6.;
    floatType C1213 =  7.;
    floatType C1222 =  8.;
    floatType C1223 =  9.;
    floatType C1233 = 10.;
    floatType C1313 = 11.;
    floatType C1322 = 12.;
    floatType C1323 = 13.;
    floatType C1333 = 14.;
    floatType C2222 = 15.;
    floatType C2223 = 16.;
    floatType C2233 = 17.;
    floatType C2323 = 18.;
    floatType C2333 = 19.;
    floatType C3333 = 20.;

    // Store the resulting stiffness tensor
    floatMatrix stiffness_tensor;

    // Fully anisotropic: 21 components
    floatVector anisotropic_parameters = { C1111, C1112, C1113, C1122, C1123, C1133, C1212, C1213, C1222, C1223, C1233,
                                           C1313, C1322, C1323, C1333, C2222, C2223, C2233, C2323, C2333, C3333 };
    floatMatrix anisotropic_answer =  {
        { C1111, C1112, C1113, C1112, C1122, C1123, C1113, C1123, C1133 },
        { C1112, C1212, C1213, C1212, C1222, C1223, C1213, C1223, C1233 },
        { C1113, C1213, C1313, C1213, C1322, C1323, C1313, C1323, C1333 },
        { C1112, C1212, C1213, C1212, C1222, C1223, C1213, C1223, C1233 },
        { C1122, C1222, C1322, C1222, C2222, C2223, C1322, C2223, C2233 },
        { C1123, C1223, C1323, C1223, C2223, C2323, C1323, C2323, C2333 },
        { C1113, C1213, C1313, C1213, C1322, C1323, C1313, C1323, C1333 },
        { C1123, C1223, C1323, C1223, C2223, C2323, C1323, C2323, C2333 },
        { C1133, C1233, C1333, C1233, C2233, C2333, C1333, C2333, C3333 }
    };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( anisotropic_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( anisotropic_answer ),
                boost::test_tools::per_element() );

    // Orthotropic: 9 components
    floatVector orthotropic_parameters = { C1111, C1122, C1133, C1212, C1313, C2222, C2233, C2323, C3333 };
    floatMatrix orthotropic_answer = {
        { C1111,    0.,    0.,    0., C1122,    0.,    0.,    0., C1133 },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        {    0.,    0., C1313,    0.,    0.,    0., C1313,    0.,    0. },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        { C1122,    0.,    0.,    0., C2222,    0.,    0.,    0., C2233 },
        {    0.,    0.,    0.,    0.,    0., C2323,    0., C2323,    0. },
        {    0.,    0., C1313,    0.,    0.,    0., C1313,    0.,    0. },
        {    0.,    0.,    0.,    0.,    0., C2323,    0., C2323,    0. },
        { C1133,    0.,    0.,    0., C2233,    0.,    0.,    0., C3333 }
    };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( orthotropic_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( orthotropic_answer ),
                boost::test_tools::per_element() );

    // Transversely isotropic/hexagonal: 5 components
    floatVector hexagonal_parameters = { C1111, C1122, C1133, C1313, C3333 };
    C1212 = 0.5 * (C1111 - C1122);  // Calculation specific to this symmetry
    floatMatrix hexagonal_answer = {
        { C1111,    0.,    0.,    0., C1122,    0.,    0.,    0., C1133 },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        {    0.,    0., C1313,    0.,    0.,    0., C1313,    0.,    0. },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        { C1122,    0.,    0.,    0., C1111,    0.,    0.,    0., C1133 },
        {    0.,    0.,    0.,    0.,    0., C1313,    0., C1313,    0. },
        {    0.,    0., C1313,    0.,    0.,    0., C1313,    0.,    0. },
        {    0.,    0.,    0.,    0.,    0., C1313,    0., C1313,    0. },
        { C1133,    0.,    0.,    0., C1133,    0.,    0.,    0., C3333 }
    };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( hexagonal_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( hexagonal_answer ),
                boost::test_tools::per_element() );
    C1212 = 6;  // Reset to match fully anisotropic index

    // Cubic symmetry: 3 components
    floatVector cubic_parameters = { C1111, C1122, C1212 };
    floatMatrix cubic_answer = {
        { C1111,    0.,    0.,    0., C1122,    0.,    0.,    0., C1122 },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        {    0.,    0., C1212,    0.,    0.,    0., C1212,    0.,    0. },
        {    0., C1212,    0., C1212,    0.,    0.,    0.,    0.,    0. },
        { C1122,    0.,    0.,    0., C1111,    0.,    0.,    0., C1122 },
        {    0.,    0.,    0.,    0.,    0., C1212,    0., C1212,    0. },
        {    0.,    0., C1212,    0.,    0.,    0., C1212,    0.,    0. },
        {    0.,    0.,    0.,    0.,    0., C1212,    0., C1212,    0. },
        { C1122,    0.,    0.,    0., C1122,    0.,    0.,    0., C1111 }
    };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( cubic_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( cubic_answer ),
                boost::test_tools::per_element() );

    floatType lambda = 12.3;
    floatType mu     = 43.4;
    floatType calc   = lambda + 2 * mu;
    floatVector isotropic_parameters = { lambda, mu };
    floatMatrix isotropic_answer =
        {
            {   calc,  0.0,  0.0,  0.0, lambda,  0.0,  0.0,  0.0, lambda },
            {    0.0,   mu,  0.0,   mu,    0.0,  0.0,  0.0,  0.0,    0.0 },
            {    0.0,  0.0,   mu,  0.0,    0.0,  0.0,   mu,  0.0,    0.0 },
            {    0.0,   mu,  0.0,   mu,    0.0,  0.0,  0.0,  0.0,    0.0 },
            { lambda,  0.0,  0.0,  0.0,   calc,  0.0,  0.0,  0.0, lambda },
            {    0.0,  0.0,  0.0,  0.0,    0.0,   mu,  0.0,   mu,    0.0 },
            {    0.0,  0.0,   mu,  0.0,    0.0,  0.0,   mu,  0.0,    0.0 },
            {    0.0,  0.0,  0.0,  0.0,    0.0,   mu,  0.0,   mu,    0.0 },
            { lambda,  0.0,  0.0,  0.0, lambda,  0.0,  0.0,  0.0,   calc }
        };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( isotropic_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( isotropic_answer ),
                boost::test_tools::per_element() );

}

BOOST_AUTO_TEST_CASE( test_evaluateEnergy ){

    unsigned int spatialDimensions = 3;

    floatVector chi = { 0.39293837, -0.42772133, -0.54629709,
                        0.10262954,  0.43893794, -0.15378708,
                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector parameters = { 12.3, 43.4 };

    floatType energy_answer = 43.983561941559444;

    floatType energy;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy ) );

    BOOST_TEST( energy == energy_answer );

    floatVector cauchyStress_answer = {-293.41192246286556,    43.439032908892734,   60.231791812665371,
                                         43.439032908892699, -224.06439311434585,    11.237268456397821,
                                         60.23179181266535,    11.237268456397864, -135.38555790789079 };


    floatVector cauchyStress;

    energy = 0;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress ) );

    BOOST_TEST( energy == energy_answer );

    BOOST_TEST( cauchyStress == cauchyStress_answer, boost::test_tools::per_element() );

    cauchyStress.clear( );

    energy = 0;

    floatVector dEnergydChi;

    floatMatrix dCauchyStressdChi;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress, dEnergydChi, dCauchyStressdChi ) );

    BOOST_TEST( energy == energy_answer );

    BOOST_TEST( cauchyStress == cauchyStress_answer, boost::test_tools::per_element() );

    // Perform tests of the gradients

    floatType eps = 1e-6;

    floatVector dEnergydChi_answer( chi.size( ), 0 );
    floatMatrix dCauchyStressdChi_answer( chi.size( ), floatVector( chi.size( ), 0 ) );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] = eps * std::abs( chi[ i ] ) + eps;

        floatType ep, em;

        floatVector cauchyStressp, cauchyStressm;

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em ) );

        dEnergydChi_answer[ i ] = ( ep - em ) / ( 2 * delta[ i ] );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep, cauchyStressp ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em, cauchyStressm ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( ( ep - em ) / ( 2 * delta[ i ] ), dEnergydChi_answer[ i ] ) );

        for ( unsigned int j = 0; j < chi.size( ); j++ ){

            dCauchyStressdChi_answer[ j ][ i ] = ( cauchyStressp[ j ] - cauchyStressm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dEnergydChi, dEnergydChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCauchyStressdChi, dCauchyStressdChi_answer ) );

    dEnergydChi = floatVector( chi.size( ), 0 );

    dCauchyStressdChi = floatMatrix( chi.size( ), floatVector( chi.size( ), 0 ) );

    floatVector d2EnergydChi2( chi.size( ) * chi.size( ), 0 );

    floatMatrix d2CauchyStressdChi2( chi.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    dEnergydChi_answer = floatVector( chi.size( ), 0 );

    dCauchyStressdChi_answer = floatMatrix( chi.size( ), floatVector( chi.size( ), 0 ) );

    floatVector d2EnergydChi2_answer( chi.size( ) * chi.size( ), 0 );

    floatMatrix d2CauchyStressdChi2_answer( chi.size( ), floatVector( chi.size( ) * chi.size( ), 0 ) );

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress, dEnergydChi, dCauchyStressdChi, d2EnergydChi2, d2CauchyStressdChi2 ) );

    BOOST_TEST( energy == energy_answer );

    BOOST_TEST( cauchyStress == cauchyStress_answer, boost::test_tools::per_element() );

    for ( unsigned int i = 0; i < chi.size( ); i++ ){

        floatVector delta( chi.size( ), 0 );

        delta[ i ] = eps * std::abs( chi[ i ] ) + eps;

        floatType ep, em;

        floatVector cauchyStressp, cauchyStressm;

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em ) );

        dEnergydChi_answer[ i ] = ( ep - em ) / ( 2 * delta[ i ] );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep, cauchyStressp ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em, cauchyStressm ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( ( ep - em ) / ( 2 * delta[ i ] ), dEnergydChi_answer[ i ] ) );

        for ( unsigned int j = 0; j < chi.size( ); j++ ){

            dCauchyStressdChi_answer[ j ][ i ] = ( cauchyStressp[ j ] - cauchyStressm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatVector dedChip, dedChim;

        floatMatrix dCauchydChip, dCauchydChim;

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep, cauchyStressp, dedChip, dCauchydChip ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em, cauchyStressm, dedChim, dCauchydChim ) );

        BOOST_CHECK( vectorTools::fuzzyEquals( ( ep - em ) / ( 2 * delta[ i ] ), dEnergydChi_answer[ i ] ) );

        for ( unsigned int j = 0; j < chi.size( ); j++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( dCauchyStressdChi_answer[ j ][ i ], ( cauchyStressp[ j ] - cauchyStressm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        for ( unsigned int k = 0; k < spatialDimensions; k++ ){
            for ( unsigned int K = 0; K < spatialDimensions; K++ ){
                d2EnergydChi2_answer[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + i ] = ( dedChip[ spatialDimensions * k + K ] - dedChim[ spatialDimensions * k + K ] ) / ( 2 * delta[ i ] );

                for ( unsigned int l = 0; l < spatialDimensions; l++ ){

                    for ( unsigned int L = 0; L < spatialDimensions; L++ ){

                        d2CauchyStressdChi2_answer[ spatialDimensions * k + K ][ spatialDimensions * spatialDimensions * spatialDimensions * l + spatialDimensions * spatialDimensions * L + i ]
                            = ( dCauchydChip[ spatialDimensions * k + K ][ spatialDimensions * l + L ] - dCauchydChim[ spatialDimensions * k + K ][ spatialDimensions * l + L ] ) / ( 2 * delta[ i ] );

                    }

                }

            }

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dEnergydChi, dEnergydChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dCauchyStressdChi, dCauchyStressdChi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2EnergydChi2, d2EnergydChi2_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( d2CauchyStressdChi2, d2CauchyStressdChi2_answer ) );

}
