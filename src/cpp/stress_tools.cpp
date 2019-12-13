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

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor with 3
         * dimensions as meanStress = frac{1}{3}*trace(\sigma)
         *
         * :param floatVector &stress: The stress tensor
         */

        //FIXME: fix mismatching argument types for trace call
        floatType trace;
        int res = vectorTools::trace(stress, trace);
        meanStress = 1./3.*trace;
         
        return NULL;
    }

    floatType calculateMeanStress(const floatVector &stress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor with 3
         * dimensions as meanStress = frac{1}{3}*trace(\sigma)
         *
         * :param floatVector &stress: The stress tensor
         */

        floatType meanStress;
        errorOut calculateMeanStressResult = calculateMeanStress(stress, meanStress);
         
        return meanStress;
    }

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain, 
                                   const floatType &previousTime, const floatVector &previousStrain, 
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         * rho^0 \Psi = 0.5*(E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ}))
         * 
         * :param const floatType &currentTime: The current time
         * :const floatVector &currentStrain: The current Green-Lagrange strain
         * :const floatType &previousTime: The previous time
         * :const floatType &currentRateModifier: The current value of the rate modifier 
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * :const floatType &previousRateModifier: The previous value of the rate modifier 
         *     which can be used for temperature effects or (potentially) other non-linear effects.
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
            factor = taui/(taui + dt*(1 - alpha)*currentRateModifier);

            //Set the indices of the previous values of the state variables
            for (unsigned int j=dim*(i-1), k=0; j<dim*i; j++, k++){
                indices[k] = j;
            }

            //Get the previous values of the state variables
            vectorTools::getValuesByIndex(previousStateVariables, indices, Xip);

            //Compute the new state-variable values
            Xic = factor*(Xip + (dt/taui)*(currentStrain*currentRateModifier + alpha*((previousStrain - Xip)*previousRateModifier - currentStrain*currentRateModifier)));

            //Add the contribution to the stress
            stress += Gi*(currentStrain - Xic);

            //Save the new value of the state variable
            currentStateVariables = vectorTools::appendVectors({currentStateVariables, Xic});
        }

        return NULL;
    }

    errorOut linearViscoelasticity(const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain, 
                                   floatVector &dstressdrateModifier){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         * rho^0 \Psi = 0.5*(E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ}))
         * 
         * :param const floatType &currentTime: The current time
         * :const floatVector &currentStrain: The current Green-Lagrange strain
         * :const floatType &previousTime: The previous time
         * :const floatVector &previousStrain: The previous value of strain
         * :const floatType &currentRateModifier: The current value of the rate modifier 
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * :const floatType &previousRateModifier: The previous value of the rate modifier 
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * :const floatVector &previousStateVariables: The previous values of the state variables
         * :const floatVector &materialParameters: The material parameters
         *     The order of the parameters is [Ginfty, taus, Gs] where
         *         Ginfty: The infinite stiffness modulus
         *         taus: The time constants
         *         Gs: The stiffness values
         * :const floatType &alpha: The integration parameter (0 for implicit, 1 for explicit)
         * :floatVector &stress: The computed stress in the reference configuration (i.e. the same configuration as the strain)
         * :floatVector &currentStateVariables: The current values of the state variables
         * :floatMatrix &dstressdstrain: The derivative of the stress w.r.t. the strain
         * :floatVector &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
         */

        errorOut lVresult = linearViscoelasticity(currentTime, currentStrain, previousTime, previousStrain, 
                                                  currentRateModifier, previousRateModifier, previousStateVariables,
                                                  materialParameters, alpha, stress, currentStateVariables);

        //Error handling
        if (lVresult){
            errorOut result = new errorNode( "linearViscoelasticity with Jacobian", "error in computation of stress");
            result->addNext(lVresult);
            return result;
        }

        //Compute the change in time
        float dt = currentTime - previousTime;

        //Compute the derivative of the stress w.r.t. the strain
        //and the rate modifier.
        dstressdrateModifier = floatVector(currentStrain.size(), 0);
        
        //Compute the "infinite" term
        floatType scalarTerm = materialParameters[0];
        floatType taui, Gi, factor;
        std::vector< unsigned int > indices(currentStrain.size(), 0);
        floatVector Xic;
        unsigned int nTerms = (materialParameters.size() - 1)/2;
        unsigned int dim = currentStrain.size();

        //Add the contributions from the other terms
        for (unsigned int i=1; i<nTerms+1; i++){
            taui = materialParameters[i];
            Gi = materialParameters[i+nTerms];

            //Set the indices of the current values of the state variables
            for (unsigned int j=dim*(i-1), k=0; j<dim*i; j++, k++){
                indices[k] = j;
            }

            //Get the previous values of the state variables
            vectorTools::getValuesByIndex(currentStateVariables, indices, Xic);

            factor = taui/(taui + dt*(1 - alpha)*currentRateModifier);

            //Add terms to the gradient w.r.t. the strain
            scalarTerm += Gi*(1 - factor*(1 - alpha)*dt/taui*currentRateModifier);

            //Add terms to the gradient w.r.t. the rate modifier
            dstressdrateModifier -= Gi*(dt*(1-alpha)/(taui + dt*(1-alpha)*currentRateModifier))*(currentStrain - Xic);
        }

        //Assemble the full gradient
        dstressdstrain = scalarTerm*vectorTools::eye<floatType>(currentStrain.size());

        return NULL;
    }

    errorOut volumetricNeoHookean(const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         * U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))
         * where J is the determinant of the deformation gradient.
         * 
         * :param const floatTyper &jacobian: The jacobian of deformation
         * :param const floatType &bulkModulus: The bulk modulus
         * :param floatType &meanStress: The meanStress NOTE: It is up to the user to determine 
         *     which configuration this mean stress is defined within. If it is the reference 
         *     configuration a mapping may be necessary going into a different configuration.
         * 
         */

        //Error handling
        if (jacobian<=0){
             return new errorNode("volumetricNeoHookean", "determinant is less than or equal zero");
        }

        //Compute the meanStress
        meanStress = 0.5*bulkModulus*(jacobian - 1/jacobian);

        return NULL;
    }

    errorOut volumetricNeoHookean(const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         * U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))
         * where J is the determinant of the deformation gradient.
         * 
         * :param const floatType &jacobian: The jacobian of deformation
         * :param const floatType &bulkModulus: The bulk modulus
         * :param floatType &meanStress: The meanStress NOTE: It is up to the user to determine 
         *     which configuration this mean stress is defined within. If it is the reference 
         *     configuration a mapping may be necessary going into a different configuration.
         * :param floatType &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian 
         *     of deformation.
         */

        //Error handling
        if (jacobian<=0){
             return new errorNode("volumetricNeoHookean", "determinant is less than or equal zero");
        }

        //Compute the meanStress
        meanStress = 0.5*bulkModulus*(jacobian - 1/jacobian);

        //Compute the derivative of the meanStress w.r.t. jacobian
        dmeanStressdJ = 0.5*bulkModulus*(1 + 1/(jacobian*jacobian));

        return NULL;
    }


    errorOut volumetricNeoHookean(const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         * U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))
         * where J is the determinant of the deformation gradient.
         * 
         * :param const floatVector &deformationGradient: The deformation gradient
         * :param const floatType &bulkModulus: The bulk modulus
         * :param floatType &meanStress: The meanStress NOTE: It is up to the user to determine 
         *     which configuration this mean stress is defined within. If it is the reference 
         *     configuration a mapping may be necessary going into a different configuration.
         * 
         */

        //Check the size of the deformation gradient
        if (deformationGradient.size() != 9){
            return new errorNode("volumetricNeoHookean", "deformation gradient must have nine terms");
        }

        //Compute the determinant of the deformation gradient
        floatType J;
        J = vectorTools::determinant(deformationGradient, 3, 3);

        return volumetricNeoHookean(J, bulkModulus, meanStress);

    }

    errorOut volumetricNeoHookean(const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         * U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))
         * where J is the determinant of the deformation gradient.
         * 
         * :param const floatVector &deformationGradient: The deformation gradient
         * :param const floatType &bulkModulus: The bulk modulus
         * :param floatType &meanStress: The meanStress NOTE: It is up to the user to determine 
         *     which configuration this mean stress is defined within. If it is the reference 
         *     configuration a mapping may be necessary going into a different configuration.
         * :param floatType &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian 
         *     of deformation.
         */

        //Check the size of the deformation gradient
        if (deformationGradient.size() != 9){
            return new errorNode("volumetricNeoHookean", "deformation gradient must have nine terms");
        }

        //Compute the determinant of the deformation gradient
        floatType J;
        J = vectorTools::determinant(deformationGradient, 3, 3);

        return volumetricNeoHookean(J, bulkModulus, meanStress, dmeanStressdJ);

    }
}
