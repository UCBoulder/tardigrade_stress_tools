//Tests for stress_tools

#include<stress_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_vector_tools
#include <boost/test/included/unit_test.hpp>

typedef constitutiveTools::errorOut errorOut;
typedef constitutiveTools::floatType floatType;
typedef constitutiveTools::floatVector floatVector;
typedef constitutiveTools::floatMatrix floatMatrix;

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

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer)
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
     *
     * :param std::ofstream &results: The output-file to write to.
     */

    //Initialize test values
    floatVector stressVector = {1., 0., 0.,
                                0., 1., 0.,
                                0., 0., 1.};
    floatMatrix stressMatrix = {{1., 0., 0.},
                                {0., 1., 0.},
                                {0., 0., 1.}};
    floatType meanStressExpected = 1.;
    floatVector jacobianVectorExpected = {1./3., 0.,    0.,
                                          0.,    1./3., 0.,
                                          0.,    0.,    1./3.};

    //Initialize test output
    floatType meanStress;
    floatVector jacobianVector(stressVector.size());

    //Test for correct mean stress calculation from stressVector with pointer output
    meanStress = 0.;
    errorOut result = stressTools::calculateMeanStress(stressVector, meanStress);
    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, meanStressExpected) );

    //Test for correct mean stress calculation from stressVector without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress(stressVector);
    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, meanStressExpected) );

    //Test for correct mean stress calculation from stressMatrix with pointer output
    meanStress = 0.;
    result = stressTools::calculateMeanStress(stressMatrix, meanStress);
    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, meanStressExpected) );

    //Test for correct mean stress calculation from stressMatrix without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress(stressMatrix);
    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, meanStressExpected) );

    //Test for correct mean stress and jacobian calculation
    meanStress = 0.;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    result = stressTools::calculateMeanStress(stressVector, meanStress, jacobianVector);
    if (!vectorTools::fuzzyEquals(meanStress, meanStressExpected) || 
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected)){
        results << "testCalculateMeanStress (test 5) & False\n";
        return 1;
    }

    results << "testCalculateMeanStress & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testCalculateDeviatoricStress ){
    /*!
     * Test the deviatoric stress calculation
     *
     * :param std::ofstream &results: The output-file to write to.
     */
    floatVector stressVector = {1., 0., 0.,
                                0., 1., 0.,
                                0., 0., 1.};
    floatVector expectedVector = {0., 0., 0.,
                                  0., 0., 0.,
                                  0., 0., 0.};
    floatVector deviatoricVector(stressVector.size(), 0.), deviatoricVectorJ;
    floatMatrix jacobian, jacobian2;
    floatType eps = 1e-6;
    errorOut result;

    //Test computation of deviatoric tensor in row major format
    std::fill(deviatoricVector.begin(), deviatoricVector.end(), 0.);
    result = stressTools::calculateDeviatoricStress(stressVector, deviatoricVector);

    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals(expectedVector, deviatoricVector) );

    std::fill(deviatoricVector.begin(), deviatoricVector.end(), 0.);

    deviatoricVector = stressTools::calculateDeviatoricStress(stressVector);

    BOOST_CHECK( vectorTools::fuzzyEquals(expectedVector, deviatoricVector) );

    //Test the computation of the jacobian
    result = stressTools::calculateDeviatoricStress(stressVector, deviatoricVectorJ, jacobian);

    BOOST_CHECK( ! result );

    BOOST_CHECK( vectorTools::fuzzyEquals(expectedVector, deviatoricVectorJ) );

    //Test the jacobian
    for (unsigned int i=0; i<stressVector.size(); i++){
        floatVector delta(stressVector.size(), 0);
        delta[i] = eps * fabs(stressVector[i]) + eps;
        
        result = stressTools::calculateDeviatoricStress(stressVector + delta, deviatoricVectorJ, jacobian);

        BOOST_CHECK( ! result );

        floatVector grad = (deviatoricVectorJ - deviatoricVector) / delta[i];

        for (unsigned int j=0; j<grad.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(grad[j], jacobian[j][i]) );
        }

    }

    deviatoricVector = stressTools::calculateDeviatoricStress(stressVector, jacobian2);

    BOOST_CHECK( vectorTools::fuzzyEquals(deviatoricVector, expectedVector) );

    BOOST_CHECK( vectorTools::fuzzyEquals(jacobian, jacobian2) );

    results << "testCalculateDeviatoricStress & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testCalculateVonMisesStress ){
    /*!
     * Test the von Mises stress calculation
     *
     * :param std::ofstream &results: The output-file to write to.
     */

    //Initialize test values
    floatVector stressVector = {1., 1., 1.,
                                1., 1., 1.,
                                1., 1., 1.};
    floatVector jacobianVectorExpected = {0.,    1./2., 1./2.,
                                          1./2., 0.,    1./2.,
                                          1./2., 1./2., 0.};
    floatType vonMisesExpected = 3.0;

    //Initialize test output
    floatType vonMises;
    floatVector jacobianVector(stressVector.size());

    //Test computation of vonMises stress from row major stress tensor
    vonMises = 0.;
    stressTools::calculateVonMisesStress(stressVector, vonMises);
    BOOST_CHECK( vectorTools::fuzzyEquals(vonMises, vonMisesExpected) );

    vonMises = 0.;
    vonMises = stressTools::calculateVonMisesStress(stressVector);
    BOOST_CHECK( vectorTools::fuzzyEquals(vonMises, vonMisesExpected) );

    //Test computation of vonMises stress and jacobian
    vonMises = 0.;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    stressTools::calculateVonMisesStress(stressVector, vonMises, jacobianVector);
    BOOST_CHECK( vectorTools::fuzzyEquals(vonMises, vonMisesExpected) &&
                 vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) );

    results << "testCalculateVonMisesStress & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testDruckerPragerSurface ){
    /*!
     * Test the Drucker-Prager yield criterion calculation.
     *
     * :param std::ofstream &results: The output-file to write to.
     */

    //Declare test input variables
    floatVector stressVector = {1., 1., 1.,
                                1., 1., 1.,
                                1., 1., 1.};
    floatType vonMises = 3.;
    floatType meanStress = 1.;
    floatType A = 1.;
    floatType B = 0.;
    floatVector dpParam = {A, B};

    floatType dpYieldExpected = 2.;
    floatVector jacobianVectorExpected = {-1./3.,  1./2.,  1./2.,
                                           1./2., -1./3.,  1./2.,
                                           1./2.,  1./2., -1./3.};
    floatVector unitDirectionVectorExpected = {-1./3.,  1./2.,  1./2.,
                                                1./2., -1./3.,  1./2.,
                                                1./2.,  1./2., -1./3.};
    unitDirectionVectorExpected /= sqrt(1.5 + 1./3);

    //Declare internal testing variables
    errorOut error;
    floatType eps;
    floatVector delta(stressVector.size(), 0.);
    floatVector gradCol(stressVector.size(), 0);
    
    //Declare test output variables
    floatType dpYield;
    floatVector jacobianVector(stressVector.size(), 0.);
    floatVector unitDirectionVector(stressVector.size(), 0.);
    floatVector jacobianVectorJ(stressVector.size(), 0.);
    floatMatrix djacobiandstress;
    floatMatrix unitDirectionJacobian;
    floatVector unitDirectionVectorJ(stressVector.size(), 0.);

    //Test computation of DP yield criterion from vonMises and meanStress
    dpYield = 0;
    stressTools::druckerPragerSurface(vonMises, meanStress, A, B, dpYield); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    stressTools::druckerPragerSurface(vonMises, meanStress, dpParam, dpYield); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(vonMises, meanStress, A, B); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(vonMises, meanStress, dpParam); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    //Test computation of DP yield criterion from row major stress tensor 
    dpYield = 0;
    stressTools::druckerPragerSurface(stressVector, A, B, dpYield); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    stressTools::druckerPragerSurface(stressVector, dpParam, dpYield); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(stressVector, A, B); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(stressVector, dpParam); 
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) );

    //Test computation of DP yield and jacobian from row major stress tensor
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    stressTools::druckerPragerSurface(stressVector, A, B, dpYield, jacobianVector);
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) &&
                 vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) );

    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    stressTools::druckerPragerSurface(stressVector, dpParam, dpYield, jacobianVector);
    BOOST_CHECK( vectorTools::fuzzyEquals(dpYield, dpYieldExpected) &&
                 vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) );

    //Test computation of DP yield, jacobian, and unit normal from row major stress tensor
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    std::fill(unitDirectionVector.begin(), unitDirectionVector.end(), 0.);
    stressTools::druckerPragerSurface(stressVector, A, B, dpYield, jacobianVector, unitDirectionVector);
    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) ||
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) ||
        !vectorTools::fuzzyEquals(unitDirectionVector, unitDirectionVectorExpected)){
        results << "testDruckerPragerSurface (test 11) & False\n";
        return 1;
    }

    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    std::fill(unitDirectionVector.begin(), unitDirectionVector.end(), 0.);
    stressTools::druckerPragerSurface(stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector);
    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) ||
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) ||
        !vectorTools::fuzzyEquals(unitDirectionVector, unitDirectionVectorExpected)){
        results << "testDruckerPragerSurface (test 12) & False\n";
        return 1;
    }

    //Test the computation of the DP yield, jacobian, and the derivative of the jacobian 
    //w.r.t. the stress
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    std::fill(jacobianVectorJ.begin(), jacobianVectorJ.end(), 0.);
    for (unsigned int i=0; i<djacobiandstress.size(); i++){
        std::fill(djacobiandstress[i].begin(), djacobiandstress[i].end(), 0.);
    }

    error = stressTools::druckerPragerSurface(stressVector, A, B, dpYield, jacobianVector, djacobiandstress);

    BOOST_CHECK( ! error );

    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) || 
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected)){
        results << "testDruckerPragerSurface (test 13) & False\n";
        return 1;
    }

    eps = 1e-6;
    for (unsigned int i=0; i<stressVector.size(); i++){
        std::fill(delta.begin(), delta.end(), 0.);
        delta[i] = eps*fabs(stressVector[i]) + eps;

        error = stressTools::druckerPragerSurface(stressVector + delta, A, B, dpYield, jacobianVectorJ);

        BOOST_CHECK( ! error );

        std::fill(gradCol.begin(), gradCol.end(), 0.);
        gradCol = (jacobianVectorJ - jacobianVector)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(djacobiandstress[j][i], gradCol[j]) );
        }
    }

    //Test the computation of the DP yield, jacobian, and the derivative of the jacobian 
    //w.r.t. the stress for the parameter vector interface
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    std::fill(jacobianVectorJ.begin(), jacobianVectorJ.end(), 0.);
    for (unsigned int i=0; i<djacobiandstress.size(); i++){
        std::fill(djacobiandstress[i].begin(), djacobiandstress[i].end(), 0.);
    }

    error = stressTools::druckerPragerSurface(stressVector, dpParam, dpYield, jacobianVector, djacobiandstress);

    BOOST_CHECK( ! error );

    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) || 
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected)){
        results << "testDruckerPragerSurface (test 15) & False\n";
        return 1;
    }

    eps = 1e-6;
    for (unsigned int i=0; i<stressVector.size(); i++){
        std::fill(delta.begin(), delta.end(), 0.);
        delta[i] = eps*fabs(stressVector[i]) + eps;

        error = stressTools::druckerPragerSurface(stressVector + delta, dpParam, dpYield, jacobianVectorJ);

        BOOST_CHECK( ! error );

        std::fill(gradCol.begin(), gradCol.end(), 0.);
        gradCol = (jacobianVectorJ - jacobianVector)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(djacobiandstress[j][i], gradCol[j]) );
        }
    }

    //Test the computation of the DP yield, jacobian, unit direction, and the jacobian 
    //of the unit direction jacobian.
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    for (unsigned int i=0; i<unitDirectionJacobian.size(); i++){
        std::fill(unitDirectionJacobian[i].begin(), unitDirectionJacobian[i].end(), 0.);
    }
    std::fill(unitDirectionVectorJ.begin(), unitDirectionVectorJ.end(), 0.);
    error = stressTools::druckerPragerSurface(stressVector, A, B, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian);

    BOOST_CHECK( ! error );

    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) || 
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) ||
        !vectorTools::fuzzyEquals(unitDirectionVector, unitDirectionVectorExpected)){
        results << "testDruckerPragerSurface (test 17) & False\n";
        return 1;
    }

    for (unsigned int i=0; i<stressVector.size(); i++){
        floatVector delta(stressVector.size(), 0);
        delta[i] = eps*fabs(stressVector[i]) + eps;

        error = stressTools::druckerPragerSurface(stressVector + delta, A, B, dpYield, jacobianVector, unitDirectionVectorJ);

        BOOST_CHECK( ! error );

        std::fill(gradCol.begin(), gradCol.end(), 0.);
        gradCol = (unitDirectionVectorJ - unitDirectionVector)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(unitDirectionJacobian[j][i], gradCol[j]) );
        }
    }

    //Test the computation of the DP yield, jacobian, unit direction, and the jacobian 
    //of the unit direction jacobian from the parameter vector interface
    dpYield = 0;
    std::fill(jacobianVector.begin(), jacobianVector.end(), 0.);
    for (unsigned int i=0; i<unitDirectionJacobian.size(); i++){
        std::fill(unitDirectionJacobian[i].begin(), unitDirectionJacobian[i].end(), 0.);
    }
    std::fill(unitDirectionVectorJ.begin(), unitDirectionVectorJ.end(), 0.);
    error = stressTools::druckerPragerSurface(stressVector, dpParam, dpYield, jacobianVector, unitDirectionVector, unitDirectionJacobian);

    BOOST_CHECK( ! error );

    if (!vectorTools::fuzzyEquals(dpYield, dpYieldExpected) || 
        !vectorTools::fuzzyEquals(jacobianVector, jacobianVectorExpected) ||
        !vectorTools::fuzzyEquals(unitDirectionVector, unitDirectionVectorExpected)){
        results << "testDruckerPragerSurface (test 19) & False\n";
        return 1;
    }

    for (unsigned int i=0; i<stressVector.size(); i++){
        floatVector delta(stressVector.size(), 0);
        delta[i] = eps*fabs(stressVector[i]) + eps;

        error = stressTools::druckerPragerSurface(stressVector + delta, dpParam, dpYield, jacobianVector, unitDirectionVectorJ);

        BOOST_CHECK( ! error );
  
        std::fill(gradCol.begin(), gradCol.end(), 0.);
        gradCol = (unitDirectionVectorJ - unitDirectionVector)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(unitDirectionJacobian[j][i], gradCol[j]) );
        }
    }
    
    results << "testDruckerPragerSurface & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testLinearViscoelasticity ){
    /*!
     * Test the implementation of linear finite-deformation
     * viscoelasticity.
     *
     * :param std::ofstream &results: The output-file to write to.
     */

    floatType previousTime = 0.6;
    floatType currentTime = 23.8;

    floatVector previousStrain = {3.03768940e-01,  4.54626930e-17, -3.71231060e-01,
                                  4.54626930e-17,  2.00000000e-01, -1.03633416e-16,
                                 -3.71231060e-01, -1.03633416e-16,  1.04623106e+00};

    floatVector currentStrain = {1.0163119 , -0.57654737,  0.33286978,
                                -0.57654737,  0.45526608, -0.22243347,
                                 0.33286978, -0.22243347,  0.19842203};

    floatVector previousStateVariables = {1, 2, 3,
                                          2, 4, 5,
                                          3, 5, 6,
                                          1.0, 0.3, 0.2,
                                          0.3, 2.0, 0.1,
                                          0.2, 0.1, 3.0};

    floatVector materialParameters = {100, 1., 10, 150, 200};

    floatType currentRateModifier = 2.7;
    floatType previousRateModifier = 3.5;

    floatType alpha = 0.5; //Trapezoidal rule

    floatVector stress;
    floatVector currentStateVariables;

    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                       previousTime, previousStrain,
                                       currentRateModifier, previousRateModifier,
                                       previousStateVariables,
                                       materialParameters, alpha,
                                       stress, currentStateVariables);

    //!Test for symmetry in the output stress
    if (!vectorTools::fuzzyEquals(stress[1], stress[3]) &&
        !vectorTools::fuzzyEquals(stress[2], stress[6]) &&
        !vectorTools::fuzzyEquals(stress[5], stress[7])){
        results << "testLinearViscoelasticity (test 1) & False\n";
        return 1;
    }

    //!Test for passing the state variables through properly
    currentTime = previousTime;

    errorOut res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                                      previousTime, previousStrain,
                                                      currentRateModifier, previousRateModifier,
                                                      previousStateVariables,
                                                      materialParameters, alpha,
                                                      stress, currentStateVariables);

    if (res){
        results << "testLinearViscoelasticity (test 2) & False\n";
        return 1;
    }

    BOOST_CHECK( vectorTools::fuzzyEquals(previousStateVariables, currentStateVariables) );

    //!Test to make sure the state variable evolution is occurring
    //!as expected for very large dt and alpha=0 (fully implicit)

    currentTime = 1e10;

    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier, previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, 0.,
                                             stress, currentStateVariables);

    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals(currentStrain*materialParameters[0], stress) );

    unsigned int dim = currentStrain.size();

    for (unsigned int i=0; i<currentStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector subv;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(currentStateVariables, indices, subv);
        BOOST_CHECK( vectorTools::fuzzyEquals(subv, currentStrain) );
    }

    //!Test for frame invariance
    floatVector Q = {-0.44956296, -0.88488713, -0.12193405,
                     -0.37866166,  0.31242661, -0.87120891,
                      0.80901699, -0.3454915 , -0.47552826};

    floatVector QT(Q.size(), 0);
    for (unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<3; j++){
            QT[3*j + i] = Q[3*i + j];
        }
    }

    //!Rotate the previous strain
    floatVector rotatedPreviousStrain(currentStrain.size(), 0);
    res = constitutiveTools::rotateMatrix(previousStrain, Q, rotatedPreviousStrain);

    BOOST_CHECK( ! res );

    //!Rotate the current strain
    floatVector rotatedCurrentStrain(currentStrain.size(), 0);
    res = constitutiveTools::rotateMatrix(currentStrain, Q, rotatedCurrentStrain);

    BOOST_CHECK( ! res );

    //!Rotate the previous state variables
    floatVector rotatedPreviousStateVariables;
    for (unsigned int i=0; i<previousStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector subv, rotatedSubv;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(previousStateVariables, indices, subv);
        res = constitutiveTools::rotateMatrix(subv, Q, rotatedSubv);
        BOOST_CHECK( ! res );
        rotatedPreviousStateVariables = vectorTools::appendVectors({rotatedPreviousStateVariables, rotatedSubv});
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
                                             rotatedStress, rotatedCurrentStateVariables);

    BOOST_CHECK( ! res );

    //Re-calculate the initial values
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier, previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, alpha,
                                             stress, currentStateVariables);

    BOOST_CHECK( ! res );

    floatVector stresspp;
    res = constitutiveTools::rotateMatrix(rotatedStress, QT, stresspp);
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals(stress, stresspp) );

    for (unsigned int i=0; i<rotatedCurrentStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector rCSV, CSV, CSVpp;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(rotatedCurrentStateVariables, indices, rCSV);
        vectorTools::getValuesByIndex(currentStateVariables, indices, CSV);

        res = constitutiveTools::rotateMatrix(rCSV, QT, CSVpp);
        BOOST_CHECK( ! res );
        BOOST_CHECK( vectorTools::fuzzyEquals(CSVpp, CSV) );
    }

    //Test the implementation of the jacobian
    floatType eps = 1e-6;
    floatVector deltaStress;

    floatMatrix jacobian;
    floatVector dstressdrateModifier;

    //Compute the jacobian
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier, previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, alpha,
                                             stress, currentStateVariables, jacobian,
                                             dstressdrateModifier);

    for (unsigned int i=0; i<currentStrain.size(); i++){
        floatVector deltaStrain(currentStrain.size(), 0);
        deltaStrain[i] = fabs(eps*currentStrain[i]);

        res = stressTools::linearViscoelasticity( currentTime,  currentStrain + deltaStrain,
                                                 previousTime, previousStrain,
                                                 currentRateModifier, previousRateModifier,
                                                 previousStateVariables,
                                                 materialParameters, alpha,
                                                 deltaStress, currentStateVariables);
        BOOST_CHECK( ! res );

        //Compare the values in the column to the jacobian's values
        for (unsigned int j=0; j<deltaStress.size(); j++){
            BOOST_CHECK( vectorTools::fuzzyEquals(jacobian[j][i], (deltaStress[j] - stress[j])/deltaStrain[i]) );
        }
    }

    floatVector deltaStrain(currentStrain.size(), 0);

    //Test the gradient w.r.t. the rate modifier
    res = stressTools::linearViscoelasticity( currentTime, currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier + fabs(eps*currentRateModifier), previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, alpha,
                                             deltaStress, currentStateVariables);

    BOOST_CHECK( ! res );

    if (! vectorTools::fuzzyEquals((deltaStress - stress)/fabs(eps*currentRateModifier), dstressdrateModifier)){
        std::cout << "dstressdrateModifier: "; vectorTools::print(dstressdrateModifier);
        std::cout << "answer:               "; vectorTools::print((deltaStress - stress)/fabs(eps*currentRateModifier));
        results << "testLinearViscoelasticity (test 6) & False\n";
        return 1;
    }

    //Test to make sure odd numbers of prony series terms can be passed in
    floatVector materialParametersOdd = { materialParameters[1], 1, 10, 100, 400, 300, 200 };
    floatVector previousStateVariablesOdd( 3 * 9, 0 );

    res = stressTools::linearViscoelasticity(  currentTime, currentStrain,
                                              previousTime, previousStrain,
                                              currentRateModifier, previousRateModifier,
                                              previousStateVariablesOdd,
                                              materialParametersOdd, alpha,
                                              stress, currentStateVariables );

    BOOST_CHECK( ! res  );
                                           

    results << "testLinearViscoelasticity & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testVolumetricNeoHookean ){
    /*!
     * Test the computation of the mean stress (-pressure) using a Neo-Hookean model
     *
     * :param std::ofstream &results: The output file
     */

    //Define the bulk modulus
    floatType bulkModulus = 150.;

    //Define the deformation gradient
    floatVector deformationGradient = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    //Run the test when there is no deformation
    floatType meanStress;

    errorOut res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, meanStress);
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, 0.) );

    //Test the meanStress computation as compared to the expected value
    deformationGradient = {0.39874077,  0.11561812, -0.75485222,
                           0.14034205,  0.15851022,  1.29640525,
                           0.26235075, -0.26051883,  0.45378251};

    floatType J = 0.25430054895115856;

    res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, meanStress);
    BOOST_CHECK( ! res );

    if (! vectorTools::fuzzyEquals(meanStress, 0.5*bulkModulus*(J - 1/J))){
        std::cout << "meanStress: " << meanStress << "\n";
        std::cout << "expected:   " << 0.5*bulkModulus*(J - 1)/J << "\n";
        results << "testVolumetricNeoHookean (test 2) & False\n";
        return 1;
    }

    //Test the meanStress computation subject to a rotation
    deformationGradient = {-0.2350804 ,  0.16410551, -1.13402371,
                            0.1296644 , -0.22975865, -1.03460443,
                           -0.4188584 , -0.16322821,  0.31618178};

    res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, meanStress);
    BOOST_CHECK( ! res );

    BOOST_CHECK( vectorTools::fuzzyEquals(meanStress, 0.5*bulkModulus*(J - 1/J)) );

    //Test the computation of the derivative of the mean stress w.r.t. J
    floatType jacobianMeanStress;
    floatType dmeanStressdJ;
    res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, jacobianMeanStress, dmeanStressdJ);
    BOOST_CHECK( ! res );

    //Make sure the mean stress is identical
    BOOST_CHECK( vectorTools::fuzzyEquals(jacobianMeanStress, meanStress) );

    //Make sure the jacobian is correct
    floatType ms;
    floatType eps = 1e-6;
    floatVector dmeanStressdF(deformationGradient.size(), 0);
    for (unsigned int i=0; i<deformationGradient.size(); i++){
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = fabs(eps*deformationGradient[i]);
        stressTools::volumetricNeoHookean(deformationGradient + delta, bulkModulus, ms);

        dmeanStressdF[i] = (ms - meanStress)/delta[i];
    }

    floatVector dJdF = vectorTools::computeDDetAdJ(deformationGradient, 3, 3);

    BOOST_CHECK( vectorTools::fuzzyEquals(dmeanStressdJ*dJdF, dmeanStressdF) );

    results << "testVolumetricNeoHookean & True\n";
    return 0;

}

BOOST_AUTO_TEST_CASE( testPeryznaModel ){
    /*!
     * Test the implementation of the Peryzna style model
     * 
     * :param std::ofstream &results: The output file
     */

    floatType f = 2.;
    floatType q = 5.42;
    floatType A = 1.4;
    floatType n = 2.4;
    floatVector parameters = {n};

    floatType p;
    errorOut error = stressTools::peryznaModel(f, q, A, n, p);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(p, A*pow((f/q), n)) );

    error = stressTools::peryznaModel(f, q, A, parameters, p);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(p, A*pow((f/q), n)) );

    floatType pJ;
    floatType dpdf, dpdq, dpdA;
    error = stressTools::peryznaModel(f, q, A, n, pJ, dpdf, dpdq, dpdA);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(p, pJ) );

    floatType eps = 1e-6;
    floatType delta = eps*fabs(f) + eps;
    error = stressTools::peryznaModel(f + delta, q, A, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdf, 1e-5, 1e-5) );
    
    delta = eps*fabs(q) + eps;
    error = stressTools::peryznaModel(f, q + delta, A, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdq, 1e-5, 1e-5) );

    delta = eps*fabs(A) + eps;
    error = stressTools::peryznaModel(f, q, A + delta, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdA, 1e-5, 1e-5) );

    f = -1;
    error = stressTools::peryznaModel(f, q, A, n, p);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(p, 0.) );

    error = stressTools::peryznaModel(f, q, A, n, pJ, dpdf, dpdq, dpdA);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(p, pJ) );

    delta = eps*fabs(f) + eps;
    error = stressTools::peryznaModel(f + delta, q, A, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdf, 1e-5, 1e-5) );
    
    delta = eps*fabs(q) + eps;
    error = stressTools::peryznaModel(f, q + delta, A, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdq, 1e-5, 1e-5) );

    delta = eps*fabs(A) + eps;
    error = stressTools::peryznaModel(f, q, A + delta, n, pJ);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals((pJ - p)/delta, dpdA, 1e-5, 1e-5) );

    floatType pJv, dpdfv, dpdqv, dpdAv;
    error = stressTools::peryznaModel(f, q, A, parameters, pJv, dpdfv, dpdqv, dpdAv);

    if (!vectorTools::fuzzyEquals(pJv, p) ||
        !vectorTools::fuzzyEquals(dpdfv, dpdf) ||
        !vectorTools::fuzzyEquals(dpdqv, dpdq) ||
        !vectorTools::fuzzyEquals(dpdAv, dpdA)){
        results << "testPeryznaModel (test 12) & False\n";
    }

    results << "testPeryznaModel & True\n";
    return 0;
}

BOOST_AUTO_TEST_CASE( testLinearHardening ){
    /*!
     * Test the linear hardening function.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector stateVariables = {1, 2, 3, 4, 5};
    floatVector linearModuli = {6, 7, 8, 9, 10};
    floatType shiftFactor = 3.7;
    floatType value;

    errorOut error = stressTools::linearHardening(stateVariables, linearModuli, shiftFactor, value);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(value, 133.7) );

    floatType valueJ;
    floatVector valueJacobian;

    error = stressTools::linearHardening(stateVariables, linearModuli, shiftFactor, valueJ, valueJacobian);

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals(value, valueJ) );

    floatType eps = 1e-6;
    for (unsigned int i=0; i<stateVariables.size(); i++){
        floatVector delta(stateVariables.size(), 0);
        delta[i] = eps*fabs(stateVariables[i]) + eps;

        error = stressTools::linearHardening(stateVariables + delta, linearModuli, shiftFactor, valueJ);
        
        BOOST_CHECK( ! error );

        BOOST_CHECK( vectorTools::fuzzyEquals((valueJ - value)/delta[i], valueJacobian[i]) );
    }

    results << "testLinearHardening & True\n";
    return 0;
}
