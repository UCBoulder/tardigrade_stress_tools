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

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
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

    floatType lamb = 12.3;
    floatType mu   = 43.4;

    floatVector parameters = { 0, lamb, mu };

    floatMatrix C_answer =
        {
            { 99.1,  0.0,  0.0,  0.0, 12.3,  0.0,  0.0,  0.0, 12.3 },
            {  0.0,  0.0,  0.0, 86.8,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 86.8,  0.0,  0.0 },
            {  0.0, 86.8,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            { 12.3,  0.0,  0.0,  0.0, 99.1,  0.0,  0.0,  0.0, 12.3 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 86.8,  0.0 },
            {  0.0,  0.0, 86.8,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0, 86.8,  0.0,  0.0,  0.0 },
            { 12.3,  0.0,  0.0,  0.0, 12.3,  0.0,  0.0,  0.0, 99.1 }
        };

    floatMatrix C;

    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( parameters, C ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( C, C_answer ) );

}

BOOST_AUTO_TEST_CASE( test_evaluateEnergy ){

    unsigned int spatialDimensions = 3;

    floatVector chi = { 0.39293837, -0.42772133, -0.54629709,
                        0.10262954,  0.43893794, -0.15378708,
                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector parameters = { 0., 12.3, 43.4 };

    floatType energy_answer = 43.98356158963631;

    floatType energy;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

    floatVector cauchyStress_answer = { -293.41192005,   43.43903202,   60.23179078,
                                          43.43903204, -224.0643919 ,   11.23726669,
                                          60.23179078,   11.23726668, -135.38556235 };

    floatVector cauchyStress;

    energy = 0;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStress, cauchyStress_answer ) );

    cauchyStress.clear( );

    energy = 0;

    floatVector dEnergydChi;

    floatMatrix dCauchyStressdChi;

    BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi, parameters, energy, cauchyStress, dEnergydChi, dCauchyStressdChi ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStress, cauchyStress_answer ) );

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

    BOOST_CHECK( vectorTools::fuzzyEquals( energy, energy_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStress, cauchyStress_answer ) );

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
