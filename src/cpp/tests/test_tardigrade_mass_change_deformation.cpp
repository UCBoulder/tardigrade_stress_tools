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

typedef tardigradeStressTools::massChangeDeformation::floatType         floatType; //!< Redefinition for the float type
typedef tardigradeStressTools::massChangeDeformation::vector3d          vector3d; //!< Redefinition for a 3d vector
typedef tardigradeStressTools::massChangeDeformation::secondOrderTensor secondOrderTensor; //!< Redefinition for the second order tensor
typedef tardigradeStressTools::massChangeDeformation::fourthOrderTensor fourthOrderTensor; //!< Redefinition for the fourth order tensor

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

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_getJA, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_getN, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    std::array< floatType, 9 > nt   = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    std::array< floatType, 9 > ntp1 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 > massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters );

    BOOST_TEST( nt   == *massChange.get_nt( ) );

    BOOST_TEST( ntp1 == *massChange.get_ntp1( ) );

}

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_getGammaRHS, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >{

        public:

            using tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >::massChangeDeformationBase;

            const floatType gammatp1 = 3.46;

            floatType JAt = 2.78;

            floatType JAtp1 = 4.13;

            secondOrderTensor nt   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            secondOrderTensor ntp1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            virtual void setJAt( ) override{
                set_JAt( JAt );
            }

            virtual void setJAtp1( ) override{
                set_JAtp1( JAtp1 );
            }

            virtual void setNt( ) override{
                set_nt( nt );
            }

            virtual void setNtp1( ) override{
                set_ntp1( ntp1 );
            }

    };

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

    floatType gammaRHS = -5.619672743296621;

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

    massChange.set_gammatp1( massChange.gammatp1 );
    BOOST_TEST( gammaRHS   == *massChange.get_gammaRHS( ) );

    floatType eps = 1e-6;

    floatType gammaLHS = 0;

    floatType dGammaRHSdJAtp1 = 0;

    secondOrderTensor dGammaRHSdNtp1;
    std::fill( std::begin( dGammaRHSdNtp1 ), std::end( dGammaRHSdNtp1 ), 0 );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( massChange.gammatp1 ) + eps;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

        floatType xp = massChange.gammatp1 + delta; 
        floatType xm = massChange.gammatp1 - delta; 

        massChangep.set_gammatp1( xp );
        massChangem.set_gammatp1( xm );

        floatType rp = *massChangep.get_gammaRHS( );
        floatType rm = *massChangem.get_gammaRHS( );

        gammaLHS = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( gammaLHS == *massChange.get_gammaLHS( ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( massChange.JAtp1 ) + eps;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_gammatp1( massChange.gammatp1 );
        massChangem.set_gammatp1( massChange.gammatp1 );

        massChangep.JAtp1 += delta;
        massChangem.JAtp1 -= delta;

        floatType rp = *massChangep.get_gammaRHS( );
        floatType rm = *massChangem.get_gammaRHS( );

        dGammaRHSdJAtp1 = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dGammaRHSdJAtp1 == *massChange.get_dGammaRHSdJAtp1( ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( massChange.ntp1[ i ] ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        vp[i] += delta;
        vm[i] -= delta;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_gammatp1( massChange.gammatp1 );
        massChangem.set_gammatp1( massChange.gammatp1 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        floatType rp = *massChangep.get_gammaRHS( );
        floatType rm = *massChangem.get_gammaRHS( );

        dGammaRHSdNtp1[ i ] = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dGammaRHSdNtp1 == *massChange.get_dGammaRHSdNtp1( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_solveGamma, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >{

        public:

            using tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >::massChangeDeformationBase;

            secondOrderTensor nt   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            secondOrderTensor ntp1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            virtual void setNt( ) override{
                set_nt( nt );
            }

            virtual void setNtp1( ) override{
                set_ntp1( ntp1 );
            }

            virtual void solveGamma_test( ){

                solveGammatp1( );

            }

    };

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

    floatType gamma_answer = -0.30717013919187797;

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

    massChange.solveGamma_test( );

    BOOST_TEST( gamma_answer == *massChange.get_gammatp1( ) );

    floatType eps = 1e-6;

    floatType dGammatp1dCtp1;

    floatType dGammatp1dRhotp1;

    secondOrderTensor dGammatp1dNtp1;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( ctp1 ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        massChangeMock massChangep( dt, At, ct, ctp1 + delta, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1 - delta, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        massChangep.solveGamma_test( );
        massChangem.solveGamma_test( );

        floatType rp = *massChangep.get_gammatp1( );
        floatType rm = *massChangem.get_gammatp1( );

        dGammatp1dCtp1 = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dGammatp1dCtp1 == *massChange.get_dGammatp1dCtp1( ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( rhotp1 ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1 + delta, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1 - delta, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        massChangep.solveGamma_test( );
        massChangem.solveGamma_test( );

        floatType rp = *massChangep.get_gammatp1( );
        floatType rm = *massChangem.get_gammatp1( );

        dGammatp1dRhotp1 = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dGammatp1dRhotp1 == *massChange.get_dGammatp1dRhotp1( ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( massChange.ntp1[ i ] ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        vp[i] += delta;
        vm[i] -= delta;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        massChangep.solveGamma_test( );
        massChangem.solveGamma_test( );

        floatType rp = *massChangep.get_gammatp1( );
        floatType rm = *massChangem.get_gammatp1( );

        dGammatp1dNtp1[ i ] = ( rp - rm ) / ( 2 * delta );

    }

    BOOST_TEST( dGammatp1dNtp1 == *massChange.get_dGammatp1dNtp1( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_computeMassDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >{

        public:

            using tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >::massChangeDeformationBase;



            secondOrderTensor nt   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            secondOrderTensor ntp1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatType convergedGammatp1 = 2.76;

            floatType dGammatp1dCtp1 = 12;
            floatType dGammatp1dRhotp1 = 12;

            secondOrderTensor dGammatp1dNtp1;

            virtual void setNt( ) override{
                set_nt( nt );
            }

            virtual void setNtp1( ) override{
                set_ntp1( ntp1 );
            }

            virtual void setConvergedGammatp1( ) override{
                set_convergedGammatp1( convergedGammatp1 );
            }

            virtual void setdGammatp1dCtp1( ) override{
                set_dGammatp1dCtp1( dGammatp1dCtp1 );
            }

            virtual void setdGammatp1dRhotp1( ) override{
                set_dGammatp1dCtp1( dGammatp1dRhotp1 );
            }

            virtual void setdGammatp1dNtp1( ) override{
                set_dGammatp1dNtp1( dGammatp1dNtp1 );
            }

    };

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

    secondOrderTensor Atp1 = { 0.66204202, -0.29969128, -0.20888693, -0.47695429,  0.50294775,
       -0.59785879, -0.53909912, -0.77874307,  0.01740997 };

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

    BOOST_TEST( Atp1 == *massChange.get_Atp1( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_massChangeDeformationBase_computeMassDeformation2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >{

        public:

            using tardigradeStressTools::massChangeDeformation::massChangeDeformationBase< 2 >::massChangeDeformationBase;



            secondOrderTensor nt   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            secondOrderTensor ntp1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            virtual void setNt( ) override{
                set_nt( nt );
            }

            virtual void setNtp1( ) override{
                set_ntp1( ntp1 );
            }

    };

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

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

    floatType eps = 1e-6;

    secondOrderTensor dAtp1dCtp1;

    secondOrderTensor dAtp1dRhotp1;

    fourthOrderTensor dAtp1dNtp1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( ctp1 ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        massChangeMock massChangep( dt, At, ct, ctp1 + delta, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1 - delta, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        secondOrderTensor rp = *massChangep.get_Atp1( );
        secondOrderTensor rm = *massChangem.get_Atp1( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dAtp1dCtp1[ j ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dAtp1dCtp1 == *massChange.get_dAtp1dCtp1( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( rhotp1 ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1 + delta, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1 - delta, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        secondOrderTensor rp = *massChangep.get_Atp1( );
        secondOrderTensor rm = *massChangem.get_Atp1( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dAtp1dRhotp1[ j ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dAtp1dRhotp1 == *massChange.get_dAtp1dRhotp1( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( massChange.ntp1[ i ] ) + eps;

        secondOrderTensor vp = massChange.ntp1;
        secondOrderTensor vm = massChange.ntp1;

        vp[ i ] += delta;
        vm[ i ] -= delta;

        massChangeMock massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );
        massChangeMock massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, parameters, 0.67 );

        massChangep.set_ntp1( vp );
        massChangem.set_ntp1( vm );

        secondOrderTensor rp = *massChangep.get_Atp1( );
        secondOrderTensor rm = *massChangem.get_Atp1( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dAtp1dNtp1[ 9 * j + i ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dAtp1dNtp1 == *massChange.get_dAtp1dNtp1( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_massChangeWeightedDirection_constructor, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    vector3d vt = { 3, 4, 5 };

    vector3d vtp1 = { 4, 5, 6 };

    std::array< floatType, 1 > parameters = { 0.7 };

    floatType alpha = 0.53;

    tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1, parameters, alpha );

    BOOST_TEST( dt         == *massChange.get_dt( ) );

    BOOST_TEST( At         == *massChange.get_At( ), CHECK_PER_ELEMENT );

    BOOST_TEST( ct         == *massChange.get_ct( ) );

    BOOST_TEST( ctp1       == *massChange.get_ctp1( ) );

    BOOST_TEST( rhot       == *massChange.get_rhot( ) );

    BOOST_TEST( rhotp1     == *massChange.get_rhotp1( ) );

    BOOST_TEST( gammat     == *massChange.get_gammat( ) );

    BOOST_TEST( parameters == *massChange.get_parameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( alpha      == *massChange.get_alpha( ) );

    BOOST_TEST( parameters[ 0 ] == *massChange.get_d( ) );

    BOOST_TEST( vt == *massChange.get_vt( ) );

    BOOST_TEST( vtp1 == *massChange.get_vtp1( ) );

}

BOOST_AUTO_TEST_CASE( test_massChangeWeightedDirection_dir, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    vector3d vt = { 3, 4, 5 };

    vector3d vtp1 = { 4, 5, 6 };

    floatType norm_vt = 7.0710678118654755;

    floatType norm_vtp1 = 8.774964387392123;

    std::array< floatType, 1 > parameters = { 0.7 };

    floatType alpha = 0.53;

    vector3d dirt = { 0.42426407, 0.56568542, 0.70710678 };

    vector3d dirtp1 = { 0.45584231, 0.56980288, 0.68376346 };

    tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1, parameters, alpha );

    BOOST_TEST( dirt == *massChange.get_dirt( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dirtp1 == *massChange.get_dirtp1( ), CHECK_PER_ELEMENT );

    BOOST_TEST( norm_vt = *massChange.get_normvt( ) );

    BOOST_TEST( norm_vtp1 = *massChange.get_normvtp1( ) );

    floatType eps = 1e-6;

    secondOrderTensor dDirtp1dVtp1;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( vtp1[ i ] ) + eps;

        vector3d vtp1p = vtp1;
        vector3d vtp1m = vtp1;

        vtp1p[ i ] += delta;
        vtp1m[ i ] -= delta;

        tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection massChangep( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1p, parameters, alpha );

        tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection massChangem( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1m, parameters, alpha );

        vector3d rp = *massChangep.get_dirtp1( );

        vector3d rm = *massChangem.get_dirtp1( );

        for ( unsigned int j = 0; j < 3; j++ ){

            dDirtp1dVtp1[ 3 * j + i ] = ( rp[ j ] - rm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dDirtp1dVtp1 == *massChange.get_dDirtp1dVtp1( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_massChangeWeightedDirection_n, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection{


        public:

            vector3d dirt   = { 1, 2, 3 };

            vector3d dirtp1 = { 2, 3, 4 };

            secondOrderTensor dDirtp1dVtp1;

            using tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection::massChangeWeightedDirection;

        public:

            virtual void setDirt( ) override{

                set_dirt( dirt );

            }

            virtual void setDirtp1( ) override{

                set_dirtp1( dirtp1 );

            }

            virtual void setdDirtp1dVtp1( ) override{

                set_dDirtp1dVtp1( dDirtp1dVtp1 );

            }

    };

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    vector3d vt = { 3, 4, 5 };

    vector3d vtp1 = { 4, 5, 6 };

    std::array< floatType, 1 > parameters = { 0.4 };

    floatType alpha = 0.53;

    secondOrderTensor nt = { 1. , 0.8, 1.2, 0.8, 2.2, 2.4, 1.2, 2.4, 4.2 };

    secondOrderTensor ntp1 = { 2.2, 2.4, 3.2, 2.4, 4.2, 4.8, 3.2, 4.8, 7. };

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1, parameters, alpha );

    BOOST_TEST( nt   == *massChange.get_nt( ), CHECK_PER_ELEMENT );

    BOOST_TEST( ntp1 == *massChange.get_ntp1( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

}

BOOST_AUTO_TEST_CASE( test_massChangeWeightedDirection_n2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection{


        public:

            vector3d dirt   = { 0, 0, 0 };

            vector3d dirtp1 = { 2, 3, 4 };

            secondOrderTensor dDirtp1dVtp1;

            using tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection::massChangeWeightedDirection;

        public:

            virtual void setNormvt( ) override{

                set_normvt( 0 );

            }

            virtual void setDirt( ) override{

                set_dirt( dirt );

            }

            virtual void setDirtp1( ) override{

                set_normvtp1( 1 );

                set_dirtp1( dirtp1 );

            }

            virtual void setdDirtp1dVtp1( ) override{

                set_dDirtp1dVtp1( dDirtp1dVtp1 );

            }

    };

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    vector3d vt = { 3, 4, 5 };

    vector3d vtp1 = { 4, 5, 6 };

    std::array< floatType, 1 > parameters = { 0.4 };

    floatType alpha = 0.53;

    secondOrderTensor nt = { 1. , 0, 0, 0, 1, 0, 0, 0, 1 };

    secondOrderTensor ntp1 = { 2.2, 2.4, 3.2, 2.4, 4.2, 4.8, 3.2, 4.8, 7. };

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1, parameters, alpha );

    BOOST_TEST( nt   == *massChange.get_nt( ), CHECK_PER_ELEMENT );

    BOOST_TEST( ntp1 == *massChange.get_ntp1( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

}

BOOST_AUTO_TEST_CASE( test_massChangeWeightedDirection_n3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class massChangeMock : public tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection{


        public:

            vector3d dirt   = { 1, 2, 3 };

            vector3d dirtp1 = { 2, 3, 4 };

            secondOrderTensor dDirtp1dVtp1;

            using tardigradeStressTools::massChangeDeformation::massChangeWeightedDirection::massChangeWeightedDirection;

        public:

            virtual void setNormvtp1( ) override{

                set_normvtp1( 0 );

            }

            virtual void setDirt( ) override{

                set_normvt( 1 );

                set_dirt( dirt );

            }

            virtual void setDirtp1( ) override{

                set_dirtp1( dirtp1 );

            }

            virtual void setdDirtp1dVtp1( ) override{

                set_dDirtp1dVtp1( dDirtp1dVtp1 );

            }

    };

    floatType dt = 1.2;

    secondOrderTensor At = { 1.03834461, -0.02177823, -0.02781574,
                             0.00522557,  1.04068676, -0.00783036,
                             0.04895802,  0.0188219 ,  1.01639564 };

    floatType ct     = 0.1;

    floatType ctp1   = 0.2;

    floatType rhot   = 1.4;

    floatType rhotp1 = 1.5;

    floatType gammat = 0.1;

    vector3d vt = { 3, 4, 5 };

    vector3d vtp1 = { 4, 5, 6 };

    std::array< floatType, 1 > parameters = { 0.4 };

    floatType alpha = 0.53;

    secondOrderTensor nt = { 1. , 0.8, 1.2, 0.8, 2.2, 2.4, 1.2, 2.4, 4.2 };

    secondOrderTensor ntp1 = { 1. , 0, 0, 0, 1, 0, 0, 0, 1 };

    massChangeMock massChange( dt, At, ct, ctp1, rhot, rhotp1, gammat, vt, vtp1, parameters, alpha );

    BOOST_TEST( nt   == *massChange.get_nt( ), CHECK_PER_ELEMENT );

    BOOST_TEST( ntp1 == *massChange.get_ntp1( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

}
