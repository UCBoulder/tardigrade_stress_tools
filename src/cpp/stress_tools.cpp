/*
===============================================================================
|                              stress_tools.cpp                               |
===============================================================================
| A collection of tools which implement and solve stress-strain relationships |
| in such a way to enable more rapid development of constitutive models which |
| have capabilities which may not be contained within a collection of         |
| of constitutive models.                                                     |
===============================================================================
*/

#include<stress_tools.h>

namespace stressTools{

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain, 
                                   const floatType &previousTime, const floatVector &previousStrain, 
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         * rho^0 \Psi = E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ})
         * 
         * :param const floatType &currentTime: The current time
         * :const floatVector &currentStrain: The current Green-Lagrange strain
         * :const floatType &previousTime: The previous time
         * :const floatVector &previousStrain: The previous value of strain
         * :const floatVector &previousStateVariables: The previous values of the state variables
         * :const floatVector &materialParameters: The material parameters
         *     The order of the parameters is [Ginfty, taus, Gs] where
         *         Ginfty: The infinite stiffness modulus
         *         taus: The time constants
         *         Gs: The stiffness values
         * :const floatType &alpha: The integration parameter (0 for implicit, 1 for explicit)
         * :floatVector &stress: The computed stress in the reference configuration (i.e. the same configuration as the strain)
         * :floatVector &currentStateVariables: The current values of the state variables
         */

        floatType dt = currentTime - previousTime;

        //Check the material parameters
        if (materialParameters.size() == 0){
            return new errorNode("linearViscoelasticity", "At least one material parameter must be provided.");
        }

        //Set the number of Prony-Series terms
        unsigned int nTerms = (materialParameters.size() - 1)/2;
        if ((nTerms % 2) != 0){
            return new errorNode("linearViscoelasticity", "An equal number of taus and Gs must be provided.");
        }

        //Set the dimension of the strain
        unsigned int dim = currentStrain.size();
        //Check the state variables
        if (previousStateVariables.size() != nTerms*dim){
            return new errorNode("linearViscoElasticity", "The number of previous state variables is not consistent with the strain size.");
        }

        //Compute the infinite stress
        stress = materialParameters[0]*currentStrain;

        //Set the initial value of the factor
        floatType factor;
        floatType taui;
        floatType Gi;
        floatVector Xip(dim, 0), Xic(dim, 0);
        currentStateVariables.resize(0);

        std::vector< unsigned int > indices(dim, 0);

        for (unsigned int i=1; i<nTerms+1; i++){
            //Get the parameters
            taui = materialParameters[i];
            Gi = materialParameters[i+nTerms];

            //Get the factor
            factor = taui/(taui + dt*(1 - alpha));

            //Set the indices of the previous values of the state variables
            for (unsigned int j=dim*i, k=0; j<dim*(i+1); j++, k++){
                indices[k] = j;
            }

            //Get the previous values of the state variables
            vectorTools::getValuesByIndex(previousStateVariables, indices, Xip);

            //Compute the new state-variable values
            Xic = factor*(Xip + (dt/taui)*(currentStrain + alpha*(previousStrain - Xip - currentStrain)));

            //Add the contribution to the stress
            stress += Gi*(currentStrain - Xic);

            //Save the new value of the state variable
            vectorTools::appendVectors({currentStateVariables, Xic});
        }

        return NULL;
    }

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         * rho^0 \Psi = E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ})
         * 
         * :param const floatType &currentTime: The current time
         * :const floatVector &currentStrain: The current Green-Lagrange strain
         * :const floatType &previousTime: The previous time
         * :const floatVector &previousStrain: The previous value of strain
         * :const floatVector &previousStateVariables: The previous values of the state variables
         * :const floatVector &materialParameters: The material parameters
         *     The order of the parameters is [Ginfty, taus, Gs] where
         *         Ginfty: The infinite stiffness modulus
         *         taus: The time constants
         *         Gs: The stiffness values
         * :const floatType &alpha: The integration parameter (0 for implicit, 1 for explicit)
         * :floatVector &stress: The computed stress in the reference configuration (i.e. the same configuration as the strain)
         * :floatVector &currentStateVariables: The current values of the state variables
         * :floatMatrix &dstressdstrain: The derivative of the stress in the reference configuration 
         *     w.r.t. the strain in the reference configuration.
         */

        errorOut lVresult = linearViscoelasticity(currentTime, currentStrain, previousTime, previousStrain, previousStateVariables,
            materialParameters, alpha, stress, currentStateVariables);

        //Error handling
        if (! lVresult){
            errorOut result( "linearViscoelasticity with Jacobian", "error in computation of stress");
            result->addNext(lVresult);
            return result;
        }

        //Compute the change in time
        float dt = currentTime - previousTime;

        //Compute the jacobian
        floatMatrix eye = vectorTools::eye<floatType>(currentStrain.size());
        
        //Compute the "infinite" term
        dstressdstrain = materialParameters[0]*eye;
        floatType taui, Gi;
        unsigned int nTerms = (materialParameters.size() - 1)/2;

        //Add the contributions from the other terms
        for (unsigned int i=1; i<nTerms+1; i++){
            taui = materialParameters[i];
            Gi = materialParameters[i+nTerms];
            dstressdstrain += Gi*(1 - dt/(taui - dt*(1-alpha))*(1-alpha))*eye;
        }

        return NULL;
    }
}
