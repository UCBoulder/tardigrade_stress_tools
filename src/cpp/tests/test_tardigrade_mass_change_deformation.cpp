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

BOOST_AUTO_TEST_CASE( test_solveMassJacobian ){

    BOOST_CHECK( true );

}
