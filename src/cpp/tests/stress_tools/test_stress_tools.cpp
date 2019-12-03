//Tests for stress_tools

#include<stress_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

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

    floatType alpha = 0.5; //Trapezoidal rule

    floatVector stress;
    floatVector currentStateVariables;

    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                       previousTime, previousStrain,
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

    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                       previousTime, previousStrain,
                                       previousStateVariables,
                                       materialParameters, alpha,
                                       stress, currentStateVariables);

    if (!vectorTools::fuzzyEquals(previousStateVariables, currentStateVariables)){
        results << "testLinearViscoelasticity (test 2) & False\n";
        return 1;
    }

    //!Test to make sure the state variable evolution is occurring 
    //!as expected for very large dt and alpha=0 (fully implicit)

    currentTime = 1e10;

    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                       previousTime, previousStrain,
                                       previousStateVariables,
                                       materialParameters, 0.,
                                       stress, currentStateVariables);

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

    //!Rotate the previous strain
    floatVector rotatedPreviousStrain(currentStrain.size(), 0);
    constitutiveTools::rotateMatrix(previousStrain, Q, rotatedPreviousStrain);

    //!Rotate the current strain
    floatVector rotatedCurrentStrain(currentStrain.size(), 0);
    constitutiveTools::rotateMatrix(currentStrain, Q, rotatedCurrentStrain);

    //!Rotate the previous state variables
    floatVector rotatedPreviousStateVariables(previousStateVariables.size(), 0);
    for (unsigned int i=0; i<previousStateVariables.size()/dim; i++){
        std::vector< unsigned int > indices(dim, 0);
        floatVector subv, rotatedSubv;

        for (unsigned int j=0; j<dim; j++){
            indices[j] = i*dim + j;
        }

        vectorTools::getValuesByIndex(previousStateVariables, indices, subv);
        constitutiveTools::rotateMatrix(subv, Q, rotatedSubv);
        rotatedPreviousStateVariables = vectorTools::appendVectors({rotatedPreviousStateVariables, rotatedSubv});
    }

    currentTime = 23.8;

    floatVector rotatedStress;
    floatVector rotatedCurrentStateVariables;

    stressTools::linearViscoelasticity( currentTime,  rotatedCurrentStrain,
                                       previousTime, rotatedPreviousStrain,
                                       rotatedPreviousStateVariables,
                                       materialParameters, alpha,
                                       rotatedStress, rotatedCurrentStateVariables);

    //Calculate the initial values
    stressTools::linearViscoelasticity( currentTime,  currentStrain,
                                       previousTime, previousStrain,
                                       previousStateVariables,
                                       materialParameters, alpha,
                                       stress, currentStateVariables);

    std::cout << "rotatedStress: "; vectorTools::print(rotatedStress);
    std::cout << "stress: "; vectorTools::print(stress);

    results << "testLinearViscoelasticity & True\n";
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
    testLinearViscoelasticity(results);

    //Close the results file
    results.close();

    return 0;
}
