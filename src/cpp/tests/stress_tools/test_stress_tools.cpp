//Tests for stress_tools

#include<stress_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

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

int testCalculateMeanStress(std::ofstream &results){
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
    errorOut result;

    //Test for correct mean stress calculation from stressVector with pointer output
    meanStress = 0.;
    result = stressTools::calculateMeanStress(stressVector, meanStress);
    if (!vectorTools::fuzzyEquals(meanStress, meanStressExpected)){
        results << "testCalculateMeanStress (test 1) & False\n";
        return 1;
    }

    //Test for correct mean stress calculation from stressVector without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress(stressVector);
    if (!vectorTools::fuzzyEquals(meanStress, meanStressExpected)){
        results << "testCalculateMeanStress (test 2) & False\n";
        return 1;
    }

    //Test for correct mean stress calculation from stressMatrix with pointer output
    meanStress = 0.;
    result = stressTools::calculateMeanStress(stressMatrix, meanStress);
    if (!vectorTools::fuzzyEquals(meanStress, meanStressExpected)){
        results << "testCalculateMeanStress (test 3) & False\n";
        return 1;
    }

    //Test for correct mean stress calculation from stressMatrix without pointer output
    meanStress = 0.;
    meanStress = stressTools::calculateMeanStress(stressMatrix);
    if (!vectorTools::fuzzyEquals(meanStress, meanStressExpected)){
        results << "testCalculateMeanStress (test 4) & False\n";
        return 1;
    }

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

int testCalculateDeviatoricStress(std::ofstream &results){
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
    floatVector deviatoricVector(stressVector.size(), 0.);
    errorOut result;

    //Test computation of deviatoric tensor in row major format
    std::fill(deviatoricVector.begin(), deviatoricVector.end(), 0.);
    result = stressTools::calculateDeviatoricStress(stressVector, deviatoricVector);
    if (!vectorTools::fuzzyEquals(expectedVector, deviatoricVector)){
        results << "testCalculateDeviatoricStress (test 1) & False\n";
        return 1;
    }

    std::fill(deviatoricVector.begin(), deviatoricVector.end(), 0.);
    deviatoricVector = stressTools::calculateDeviatoricStress(stressVector);
    if (!vectorTools::fuzzyEquals(expectedVector, deviatoricVector)){
        results << "testCalculateDeviatoricStress (test 2) & False\n";
        return 1;
    }

    results << "testCalculateDeviatoricStress & True\n";
    return 0;
}

int testCalculateVonMisesStress(std::ofstream &results){
    /*!
     * Test the von Mises stress calculation
     *
     * :param std::ofstream &results: The output-file to write to.
     */
    floatVector stressVector = {1., 1., 1.,
                                1., 1., 1.,
                                1., 1., 1.};
    floatType expected = 3.0;
    floatType vonMises;

    //Test computation of vonMises stress from row major stress tensor
    vonMises = 0.;
    stressTools::calculateVonMisesStress(stressVector, vonMises);
    if (!vectorTools::fuzzyEquals(vonMises, expected)){
        results << "testCalculateVonMisesStress (test 1) & False\n";
        return 1;
    }

    vonMises = 0.;
    vonMises = stressTools::calculateVonMisesStress(stressVector);
    if (!vectorTools::fuzzyEquals(vonMises, expected)){
        results << "testCalculateVonMisesStress (test 2) & False\n";
        return 1;
    }

    results << "testCalculateVonMisesStress & True\n";
    return 0;
}

int testDruckerPragerSurface(std::ofstream &results){
    /*!
     * Test the Drucker-Prager yield criterion calculation.
     *
     * :param std::ofstream &results: The output-file to write to.
     */

    floatVector stressVector = {1., 1., 1.,
                                1., 1., 1.,
                                1., 1., 1.};
    floatType vonMises = 3.;
    floatType meanStress = 1.;
    floatType A = 1.;

    floatType expected = 2.;
    floatType dpYield;

    //Test computation of DP yield criterion from vonMises and meanStress
    dpYield = 0;
    stressTools::druckerPragerSurface(vonMises, meanStress, A, dpYield); 
    if (!vectorTools::fuzzyEquals(dpYield, expected)){
        results << "testDruckerPragerSurface (test 1) & False\n";
        return 1;
    }

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(vonMises, meanStress, A); 
    if (!vectorTools::fuzzyEquals(dpYield, expected)){
        results << "testDruckerPragerSurface (test 2) & False\n";
        return 1;
    }

    //Test computation of DP yield criterion from row major stress tensor 
    dpYield = 0;
    stressTools::druckerPragerSurface(stressVector, A, dpYield); 
    if (!vectorTools::fuzzyEquals(dpYield, expected)){
        results << "testDruckerPragerSurface (test 3) & False\n";
        return 1;
    }

    dpYield = 0;
    dpYield = stressTools::druckerPragerSurface(stressVector, A); 
    if (!vectorTools::fuzzyEquals(dpYield, expected)){
        results << "testDruckerPragerSurface (test 4) & False\n";
        return 1;
    }

    results << "testDruckerPragerSurface & True\n";
    return 0;
}

int testLinearViscoelasticity(std::ofstream &results){
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

    if (!vectorTools::fuzzyEquals(previousStateVariables, currentStateVariables)){
        results << "testLinearViscoelasticity (test 2) & False\n";
        return 1;
    }

    //!Test to make sure the state variable evolution is occurring
    //!as expected for very large dt and alpha=0 (fully implicit)

    currentTime = 1e10;

    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier, previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, 0.,
                                             stress, currentStateVariables);

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 3) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(currentStrain*materialParameters[0], stress)){
        results << "testLinearViscoelasticity (test 3) & False\n";
        return 1;
    }

    unsigned int dim = currentStrain.size();

    for (unsigned int i=0; i<currentStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector subv;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(currentStateVariables, indices, subv);
        if (!vectorTools::fuzzyEquals(subv, currentStrain)){
            results << "testLinearViscoelasticity (test 4) & False\n";
            return 1;
        }
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

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

    //!Rotate the current strain
    floatVector rotatedCurrentStrain(currentStrain.size(), 0);
    res = constitutiveTools::rotateMatrix(currentStrain, Q, rotatedCurrentStrain);

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

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
        if (res){
            res->print();
            results << "testLinearViscoelasticity (test 5) & False\n";
            return 1;
        }
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

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

    //Re-calculate the initial values
    res = stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                             previousTime, previousStrain,
                                             currentRateModifier, previousRateModifier,
                                             previousStateVariables,
                                             materialParameters, alpha,
                                             stress, currentStateVariables);

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

    floatVector stresspp;
    res = constitutiveTools::rotateMatrix(rotatedStress, QT, stresspp);
    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(stress, stresspp)){
        results << "testLinearViscoelasticity (test 5) & False\n";
        return 1;
    }

    for (unsigned int i=0; i<rotatedCurrentStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector rCSV, CSV, CSVpp;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(rotatedCurrentStateVariables, indices, rCSV);
        vectorTools::getValuesByIndex(currentStateVariables, indices, CSV);

        res = constitutiveTools::rotateMatrix(rCSV, QT, CSVpp);
        if (res){
            res->print();
            results << "testLinearViscoelasticity (test 5) & False\n";
            return 1;
        }
        if (!vectorTools::fuzzyEquals(CSVpp, CSV)){
            results << "testLinearViscoelasticity (test 5) & False\n";
            return 1;
        }
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
        if (res){
            res->print();
            results << "testLinearViscoelasticity (test 5) & False\n";
            return 1;
        }

        //Compare the values in the column to the jacobian's values
        for (unsigned int j=0; j<deltaStress.size(); j++){
            if (! vectorTools::fuzzyEquals(jacobian[j][i], (deltaStress[j] - stress[j])/deltaStrain[i])){
                results << "testLinearViscoelasticity (test 5) & False\n";
                return 1;
            }
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

    if (res){
        res->print();
        results << "testLinearViscoelasticity (test 6) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals((deltaStress - stress)/fabs(eps*currentRateModifier), dstressdrateModifier)){
        std::cout << "dstressdrateModifier: "; vectorTools::print(dstressdrateModifier);
        std::cout << "answer:               "; vectorTools::print((deltaStress - stress)/fabs(eps*currentRateModifier));
        results << "testLinearViscoelasticity (test 6) & False\n";
        return 1;
    }

    results << "testLinearViscoelasticity & True\n";
    return 0;
}

int testVolumetricNeoHookean(std::ofstream &results){
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
    if (res){
        res->print();
        results << "testVolumetricNeoHookean (test 1) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals(meanStress, 0.)){
        results << "testVolumetricNeoHookean (test 1) & False\n";
        return 1;
    }

    //Test the meanStress computation as compared to the expected value
    deformationGradient = {0.39874077,  0.11561812, -0.75485222,
                           0.14034205,  0.15851022,  1.29640525,
                           0.26235075, -0.26051883,  0.45378251};

    floatType J = 0.25430054895115856;

    res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, meanStress);
    if (res){
        res->print();
        results << "testVolumetricNeoHookean (test 2) & False\n";
        return 1;
    }

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
    if (res){
        res->print();
        results << "testVolumetricNeoHookean (test 3) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals(meanStress, 0.5*bulkModulus*(J - 1/J))){
        results << "testVolumetricNeoHookean (test 3) & False\n";
        return 1;
    }

    //Test the computation of the derivative of the mean stress w.r.t. J
    floatType jacobianMeanStress;
    floatType dmeanStressdJ;
    res = stressTools::volumetricNeoHookean(deformationGradient, bulkModulus, jacobianMeanStress, dmeanStressdJ);
    if (res){
        res->print();
        results << "testVolumetricNeoHookean (test 4) & False\n";
        return 1;
    }

    //Make sure the mean stress is identical
    if (!vectorTools::fuzzyEquals(jacobianMeanStress, meanStress)){
        results << "testVolumetricNeoHookean (test 4) & False\n";
        return 1;
    }

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

    if (!vectorTools::fuzzyEquals(dmeanStressdJ*dJdF, dmeanStressdF)){
        results << "testVolumentricNeoHookean (test 4) & False\n";
        return 1;
    }

    results << "testVolumetricNeoHookean & True\n";
    return 0;

}

int main(){
    /*!
    The main loop which runs the tests defined in the
    accompanying functions. Each function should output
    the function name followed by & followed by True or False
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    testCalculateMeanStress(results);
    testCalculateVonMisesStress(results);
    testDruckerPragerSurface(results);
    testCalculateDeviatoricStress(results);
    testLinearViscoelasticity(results);
    testVolumetricNeoHookean(results);

    //Close the results file
    results.close();

    return 0;
}
