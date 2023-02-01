/**
  * \file test_linear_elasticity.cpp
  *
  * Tests for the linear elasticity support module
  */

#include<stress_tools.h>
#include<sstream>
#include<fstream>
#include<math.h>

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

    // Full 9x9 as a row-major vector
    floatVector eightyone_parameters = {  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,
                                         10., 11., 12., 13., 14., 15., 16., 17., 18.,
                                         19., 20., 21., 22., 23., 24., 25., 26., 27.,
                                         28., 29., 30., 31., 32., 33., 34., 35., 36.,
                                         37., 38., 39., 40., 41., 42., 43., 44., 45.,
                                         46., 47., 48., 49., 50., 51., 52., 53., 54.,
                                         55., 56., 57., 58., 59., 60., 61., 62., 63.,
                                         64., 65., 66., 67., 68., 69., 70., 71., 72.,
                                         73., 74., 75., 76., 77., 78., 79., 80., 81. };
    floatMatrix eightyone_answer = { {  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9. },
                                     { 10., 11., 12., 13., 14., 15., 16., 17., 18. },
                                     { 19., 20., 21., 22., 23., 24., 25., 26., 27. },
                                     { 28., 29., 30., 31., 32., 33., 34., 35., 36. },
                                     { 37., 38., 39., 40., 41., 42., 43., 44., 45. },
                                     { 46., 47., 48., 49., 50., 51., 52., 53., 54. },
                                     { 55., 56., 57., 58., 59., 60., 61., 62., 63. },
                                     { 64., 65., 66., 67., 68., 69., 70., 71., 72. },
                                     { 73., 74., 75., 76., 77., 78., 79., 80., 81. } };
    BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( eightyone_parameters, stiffness_tensor ) );
    BOOST_TEST( vectorTools::appendVectors( stiffness_tensor ) == vectorTools::appendVectors( eightyone_answer ),
                boost::test_tools::per_element() );

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

BOOST_AUTO_TEST_CASE( test_rotations_formReferenceStiffnessTensor, * boost::unit_test::tolerance( 1.0e-13 ) ){

    unsigned int spatialDimensions = 3;
    floatMatrix stiffnessTensor;
    floatMatrix directionCosines;
    floatVector parameters = {  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,
                               10., 11., 12., 13., 14., 15., 16., 17., 18.,
                               19., 20., 21., 22., 23., 24., 25., 26., 27.,
                               28., 29., 30., 31., 32., 33., 34., 35., 36.,
                               37., 38., 39., 40., 41., 42., 43., 44., 45.,
                               46., 47., 48., 49., 50., 51., 52., 53., 54.,
                               55., 56., 57., 58., 59., 60., 61., 62., 63.,
                               64., 65., 66., 67., 68., 69., 70., 71., 72.,
                               73., 74., 75., 76., 77., 78., 79., 80., 81. };
    floatMatrix bungeEulerAngles = {
        {   0.,     0.,   0. },
        { M_PI,     0.,   0. },
        {   0.,     0., M_PI },
        {   0.,   M_PI,   0. },
        { M_PI, M_PI_2,   0. },
        {   0., M_PI_2, M_PI }
    };
    std::vector< std::vector< std::vector< double > > > expected_stiffnessTensor = {
        { {  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9. },
          { 10., 11., 12., 13., 14., 15., 16., 17., 18. },
          { 19., 20., 21., 22., 23., 24., 25., 26., 27. },
          { 28., 29., 30., 31., 32., 33., 34., 35., 36. },
          { 37., 38., 39., 40., 41., 42., 43., 44., 45. },
          { 46., 47., 48., 49., 50., 51., 52., 53., 54. },
          { 55., 56., 57., 58., 59., 60., 61., 62., 63. },
          { 64., 65., 66., 67., 68., 69., 70., 71., 72. },
          { 73., 74., 75., 76., 77., 78., 79., 80., 81. } },
        { {   1.,   2.,  -3.,   4.,   5.,  -6.,  -7.,  -8.,   9. },
          {  10.,  11., -12.,  13.,  14., -15., -16., -17.,  18. },
          { -19., -20.,  21., -22., -23.,  24.,  25.,  26., -27. },
          {  28.,  29., -30.,  31.,  32., -33., -34., -35.,  36. },
          {  37.,  38., -39.,  40.,  41., -42., -43., -44.,  45. },
          { -46., -47.,  48., -49., -50.,  51.,  52.,  53., -54. },
          { -55., -56.,  57., -58., -59.,  60.,  61.,  62., -63. },
          { -64., -65.,  66., -67., -68.,  69.,  70.,  71., -72. },
          {  73.,  74., -75.,  76.,  77., -78., -79., -80.,  81. } },
        { {   1.,   2.,  -3.,   4.,   5.,  -6.,  -7.,  -8.,   9. },
          {  10.,  11., -12.,  13.,  14., -15., -16., -17.,  18. },
          { -19., -20.,  21., -22., -23.,  24.,  25.,  26., -27. },
          {  28.,  29., -30.,  31.,  32., -33., -34., -35.,  36. },
          {  37.,  38., -39.,  40.,  41., -42., -43., -44.,  45. },
          { -46., -47.,  48., -49., -50.,  51.,  52.,  53., -54. },
          { -55., -56.,  57., -58., -59.,  60.,  61.,  62., -63. },
          { -64., -65.,  66., -67., -68.,  69.,  70.,  71., -72. },
          {  73.,  74., -75.,  76.,  77., -78., -79., -80.,  81. } },
        { {   1.,  -2.,  -3.,  -4.,   5.,   6.,  -7.,   8.,   9. },
          { -10.,  11.,  12.,  13., -14., -15.,  16., -17., -18. },
          { -19.,  20.,  21.,  22., -23., -24.,  25., -26., -27. },
          { -28.,  29.,  30.,  31., -32., -33.,  34., -35., -36. },
          {  37., -38., -39., -40.,  41.,  42., -43.,  44.,  45. },
          {  46., -47., -48., -49.,  50.,  51., -52.,  53.,  54. },
          { -55.,  56.,  57.,  58., -59., -60.,  61., -62., -63. },
          {  64., -65., -66., -67.,  68.,  69., -70.,  71.,  72. },
          {  73., -74., -75., -76.,  77.,  78., -79.,  80.,  81. } },
        { {   1.,  -3.,  -2.,  -7.,   9.,   8.,  -4.,   6.,   5. },
          { -19.,  21.,  20.,  25., -27., -26.,  22., -24., -23. },
          { -10.,  12.,  11.,  16., -18., -17.,  13., -15., -14. },
          { -55.,  57.,  56.,  61., -63., -62.,  58., -60., -59. },
          {  73., -75., -74., -79.,  81.,  80., -76.,  78.,  77. },
          {  64., -66., -65., -70.,  72.,  71., -67.,  69.,  68. },
          { -28.,  30.,  29.,  34., -36., -35.,  31., -33., -32. },
          {  46., -48., -47., -52.,  54.,  53., -49.,  51.,  50. },
          {  37., -39., -38., -43.,  45.,  44., -40.,  42.,  41. } },
        { {  1.,  3.,  2.,  7.,  9.,  8.,  4.,  6.,  5. },
          { 19., 21., 20., 25., 27., 26., 22., 24., 23. },
          { 10., 12., 11., 16., 18., 17., 13., 15., 14. },
          { 55., 57., 56., 61., 63., 62., 58., 60., 59. },
          { 73., 75., 74., 79., 81., 80., 76., 78., 77. },
          { 64., 66., 65., 70., 72., 71., 67., 69., 68. },
          { 28., 30., 29., 34., 36., 35., 31., 33., 32. },
          { 46., 48., 47., 52., 54., 53., 49., 51., 50. },
          { 37., 39., 38., 43., 45., 44., 40., 42., 41. } },
    };


    for ( unsigned int i=0; i<bungeEulerAngles.size( ); i++ ){

        directionCosines = floatMatrix( spatialDimensions, floatVector( spatialDimensions, 0 ) );
        vectorTools::rotationMatrix( bungeEulerAngles[ i ], directionCosines );

        //Test directionCosines interface
        stiffnessTensor = floatMatrix( spatialDimensions * spatialDimensions,
                                       floatVector( spatialDimensions * spatialDimensions, 0 ) );
        BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( directionCosines, parameters,
                                                                                   stiffnessTensor ) );
        BOOST_TEST( vectorTools::appendVectors( stiffnessTensor ) ==
                        vectorTools::appendVectors( expected_stiffnessTensor[ i ] ),
                    boost::test_tools::per_element() );

        //Test bungeEulerAngles interface
        stiffnessTensor = floatMatrix( spatialDimensions * spatialDimensions,
                                       floatVector( spatialDimensions * spatialDimensions, 0 ) );
        BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( bungeEulerAngles[ i ], parameters,
                                                                                   stiffnessTensor ) );
        BOOST_TEST( vectorTools::appendVectors( stiffnessTensor ) ==
                        vectorTools::appendVectors( expected_stiffnessTensor[ i ] ),
                    boost::test_tools::per_element() );
    }

    //Test a rotation that doesn't just rearrange the stiffness tensor components
    parameters = floatVector( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 1 );
    bungeEulerAngles = {
        { M_PI_4, M_PI_4,   0. }
    };
    expected_stiffnessTensor = {
        { {0.25, 0.25, 0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  1. },
          {0.25, 0.25, 0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  1. },
          {0.5,  0.5,  1.,   0.5,  0.5,  1.,   1.,   1.,   2. },
          {0.25, 0.25, 0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  1. },
          {0.25, 0.25, 0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  1. },
          {0.5,  0.5,  1.,   0.5,  0.5,  1.,   1.,   1.,   2. },
          {0.5,  0.5,  1.,   0.5,  0.5,  1.,   1.,   1.,   2. },
          {0.5,  0.5,  1.,   0.5,  0.5,  1.,   1.,   1.,   2. },
          {1.,   1.,   2.,   1.,   1.,   2.,   2.,   2.,   4. } }
    };
    for ( unsigned int i=0; i<bungeEulerAngles.size( ); i++ ){

        directionCosines = floatMatrix( spatialDimensions, floatVector( spatialDimensions, 0 ) );
        vectorTools::rotationMatrix( bungeEulerAngles[ i ], directionCosines );

        //Test directionCosines interface
        stiffnessTensor = floatMatrix( spatialDimensions * spatialDimensions,
                                       floatVector( spatialDimensions * spatialDimensions, 0 ) );
        BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( directionCosines, parameters,
                                                                                   stiffnessTensor ) );
        BOOST_TEST( vectorTools::appendVectors( stiffnessTensor ) ==
                        vectorTools::appendVectors( expected_stiffnessTensor[ i ] ),
                    boost::test_tools::per_element() );

        //Test bungeEulerAngles interface
        stiffnessTensor = floatMatrix( spatialDimensions * spatialDimensions,
                                       floatVector( spatialDimensions * spatialDimensions, 0 ) );
        BOOST_CHECK( !stressTools::linearElasticity::formReferenceStiffnessTensor( bungeEulerAngles[ i ], parameters,
                                                                                   stiffnessTensor ) );
        BOOST_TEST( vectorTools::appendVectors( stiffnessTensor ) ==
                        vectorTools::appendVectors( expected_stiffnessTensor[ i ] ),
                    boost::test_tools::per_element() );
    }

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

        BOOST_TEST( ( ep - em ) / ( 2 * delta[ i ] ) == dEnergydChi_answer[ i ] );

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

        BOOST_TEST( ( ep - em ) / ( 2 * delta[ i ] ) == dEnergydChi_answer[ i ] );

        for ( unsigned int j = 0; j < chi.size( ); j++ ){

            dCauchyStressdChi_answer[ j ][ i ] = ( cauchyStressp[ j ] - cauchyStressm[ j ] ) / ( 2 * delta[ i ] );

        }

        floatVector dedChip, dedChim;

        floatMatrix dCauchydChip, dCauchydChim;

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi + delta, parameters, ep, cauchyStressp, dedChip, dCauchydChip ) );

        BOOST_CHECK( !stressTools::linearElasticity::evaluateEnergy( chi - delta, parameters, em, cauchyStressm, dedChim, dCauchydChim ) );

        BOOST_TEST( ( ep - em ) / ( 2 * delta[ i ] ) == dEnergydChi_answer[ i ] );

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
