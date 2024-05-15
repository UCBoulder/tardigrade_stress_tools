/**
  * \file test_mass_change_deformation.cpp
  *
  * Tests for the mass change deformation header only library
  */

#include<tardigrade_mass_change_deformation.h>
#include<sstream>
#include<fstream>
#include<math.h>

#define BOOST_TEST_MODULE test_linear_elasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeStressTools::massChangeDeformation::floatType floatType; //!< Redefinition for the float type
typedef tardigradeStressTools::massChangeDeformation::secondOrderTensor secondOrderTensor; //!< Redefinition for the second order tensor

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_constructor, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    std::array< floatType, 2 > parameters = { 1, 2 };

    floatType alpha = 0.53;

    tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, alpha );

    BOOST_TEST( dt         == *massChange.get_dt( ) );

    BOOST_TEST( At         == *massChange.get_At( ), CHECK_PER_ELEMENT );

    BOOST_TEST( ct         == *massChange.get_ct( ) );

    BOOST_TEST( ctp1       == *massChange.get_ctp1( ) );

    BOOST_TEST( rhot       == *massChange.get_rhot( ) );

    BOOST_TEST( rhotp1     == *massChange.get_rhotp1( ) );

    BOOST_TEST( gammat     == *massChange.get_gammat( ) );

    BOOST_TEST( parameters == *massChange.get_parameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( alpha      == *massChange.get_alpha( ) );

}

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_solveMassJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType JAt = 1.1;

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    std::array< floatType, 2 > parameters = { 1, 2 };

    floatType JAtp1 = 1.2468944099378887;

    tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters );

    BOOST_TEST( JAt   == *massChange.get_JAt( ) );

    BOOST_TEST( JAtp1 == *massChange.get_JAtp1( ) );

    floatType eps = 1e-3;

    floatType dJAtp1dCtp1 = 0;

    floatType dJAtp1dRhotp1 = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( ctp1 ) + eps;

        tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChangep( dt, At, ct, ctp1 + delta, rhot, rhotp1, gammat, parameters );
        tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChangem( dt, At, ct, ctp1 - delta, rhot, rhotp1, gammat, parameters );
        
        floatType rp = *massChangep.get_JAtp1( );
        floatType rm = *massChangem.get_JAtp1( );

        dJAtp1dCtp1 = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dJAtp1dCtp1 == *massChange.get_dJAtp1dCtp1( ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( ctp1 ) + eps;

        tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChangep( dt, At, ct, ctp1, rhot, rhotp1 + delta, gammat, parameters );
        tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChangem( dt, At, ct, ctp1, rhot, rhotp1 - delta, gammat, parameters );
        
        floatType rp = *massChangep.get_JAtp1( );
        floatType rm = *massChangem.get_JAtp1( );

        dJAtp1dRhotp1 = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dJAtp1dRhotp1 == *massChange.get_dJAtp1dRhotp1( ) );

}
