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
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         * meanStress = \frac{1}{3}*trace(\sigma)
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &meanStress: The mean stress scalar 
         */

        floatType trace;
        int result = vectorTools::trace(stress, trace);
        meanStress = 1./3.*trace;
         
        return NULL;
    }

    floatType calculateMeanStress(const floatVector &stress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         * meanStress = \frac{1}{3}*trace(\sigma)
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &meanStress: The mean stress scalar 
         */

        floatType meanStress;
        errorOut result = calculateMeanStress(stress, meanStress);
         
        return meanStress;
    }

    errorOut calculateMeanStress(const floatMatrix &stress, floatType &meanStress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         * meanStress = \frac{1}{3}*trace(\sigma)
         *
         * :param floatMatrix &stress: The stress tensor
         * :param floatType &meanStress: The mean stress scalar 
         */

        floatType trace;
        int result = vectorTools::trace(stress, trace);
        meanStress = 1./3.*trace;
         
        return NULL;
    }

    floatType calculateMeanStress(const floatMatrix &stress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         * meanStress = \frac{1}{3}*trace(\sigma)
         *
         * :param floatMatrix &stress: The stress tensor
         * :param floatType &meanStress: The mean stress scalar 
         */

        floatType meanStress;
        errorOut result = calculateMeanStress(stress, meanStress);
         
        return meanStress;
    }

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress, floatVector &jacobian){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         * meanStress = \frac{1}{3}*trace(\sigma)
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &meanStress: The scalar mean stress  
         * :param floatVector &jacobian: The row major jacobian w.r.t. the stress tensor
         */

        //Check vector lengths
        unsigned int length = stress.size();
        if (length != jacobian.size()){
            return new errorNode("calculateMeanStress", "The stress tensor and jacobian tensor sizes must match.");
        }

        //Calculate the mean stress
        meanStress = 0.;
        std::fill(jacobian.begin(), jacobian.end(), 0.);
        calculateMeanStress(stress, meanStress);

        //Initialize the identity matrix
        floatVector I(length, 0.);
        vectorTools::eye<floatType>(I);

        //Calculate the jacobian
        jacobian = 1./3.*I;
         
        return NULL;
    }

    errorOut calculateDeviatoricStress(const floatVector &stress, floatVector &deviatoric){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         * \sigma^{deviatoric} = \sigma - \sigma^{mean}I
         *
         * :param floatVector &stress: The stress tensor in row major format
         * :param floatVector &deviatoric: The deviatoric stress tensor in row major format
         */

        //Check vector lengths
        unsigned int length = stress.size();
        if (length != deviatoric.size()){
            return new errorNode("calculateDeviatoric", "The tensor and deviatoric tensor sizes must match.");
        }

        //Initialize the identity matrix
        floatVector I(length, 0.);
        vectorTools::eye<floatType>(I);

        //Calculate deviatoric stress tensor
        floatType meanStress = calculateMeanStress(stress);
        deviatoric = stress - meanStress*I;
    
        return NULL;
    }

    floatVector calculateDeviatoricStress(const floatVector &stress){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         * \sigma^{deviatoric} = \sigma - \sigma^{mean}I
         *
         * :param floatVector &stress: The stress tensor in row major format
         * :returns: The deviatoric stress tensor in row major format
         * :rtype: floatVector &deviatoric
         */

        floatVector deviatoric(stress.size());
        errorOut result;
        result = calculateDeviatoricStress(stress, deviatoric); 
    
        return deviatoric;
    }

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         * \sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}
         * \sigma^{deviatoric} = \sigma - \sigma^{mean}I
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &meanStress: The scalar mean stress 
         */

        floatVector deviatoric = calculateDeviatoricStress(stress);
        vonMises = std::sqrt(3./2.*vectorTools::inner(deviatoric, deviatoric));

        return NULL;
    }

    floatType calculateVonMisesStress(const floatVector &stress){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         * \sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}
         * \sigma^{deviatoric} = \sigma - \sigma^{mean}I
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &meanStress: The scalar mean stress 
         */

        floatType vonMises = 0.;
        errorOut result = calculateVonMisesStress(stress, vonMises);

        return vonMises;
    }

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises, floatVector &jacobian){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         * \sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}
         * \sigma^{deviatoric} = \sigma - \sigma^{mean}I
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &meanStress: The scalar mean stress 
         * :param floatVector &jacobian: The row major mean stress jacobian tensor w.r.t. the stress tensor
         */

        //Check vector lengths
        unsigned int length = stress.size();
        if (length != jacobian.size()){
            return new errorNode("calculateVonMisesStress", "The stress tensor and jacobian tensor sizes must match.");
        }

        //Calculate the vonMises stress
        calculateVonMisesStress(stress, vonMises);

        //Calculate the deviatoric stress
        floatVector deviatoric(stress.size(), 0.);
        calculateDeviatoricStress(stress, deviatoric);

        //Calculate the jacobian
        jacobian = 3./(2.*vonMises) * deviatoric;

        return NULL;
    }

    errorOut druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress 
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        dpYield = vonMises - A*meanStress - B;
    
        return NULL;
    }

    errorOut druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress 
         * f = \sigma^{vonMises} - dpParam[0]*\sigma^{mean} - dpParam[1]
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatVector &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        //Check Drucker-Prager parameter vector length 
        if (dpParam.size() != 2){
            return new errorNode("druckerPragerSurface", "Two parameters are required for the Drucker-Prager surface.");
        }

        druckerPragerSurface(vonMises, meanStress, dpParam[0], dpParam[1], dpYield);
    
        return NULL;
    }

    floatType druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress 
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B 
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :returns:  The Drucker-Prager yield stress/criterion/surface
         * :rtype: floatType &dpYield
         */

        floatType dpYield = 0;

        //Calculate DP yield criterion
        druckerPragerSurface(vonMises, meanStress, A, B, dpYield);
    
        return dpYield;
    }

    floatType druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress 
         * f = \sigma^{vonMises} - dpParam[0]*\sigma^{mean} - dpParam[1]
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatVector &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :returns:  The Drucker-Prager yield stress/criterion/surface
         * :rtype: floatType &dpYield
         */

        floatType dpYield = 0;
        druckerPragerSurface(vonMises, meanStress, dpParam, dpYield);
    
        return dpYield;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        //Calculate von Mises and mean stresses
        floatType vonMises, meanStress;
        vonMises = meanStress = 0.;
        calculateVonMisesStress(stress, vonMises);
        calculateMeanStress(stress, meanStress);

        //Calculate DP yield criterion
        druckerPragerSurface(vonMises, meanStress, A, B, dpYield);
    
        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatVector &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        //Check Drucker-Prager parameter vector length 
        if (dpParam.size() != 2){
            return new errorNode("druckerPragerSurface", "Two parameters are required for the Drucker-Prager surface.");
        }

        //Calculate DP yield criterion
        druckerPragerSurface(stress, dpParam[0], dpParam[1], dpYield);
    
        return NULL;
    }

    floatType druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The row major stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :returns: The Drucker-Prager yield stress/criterion/surface
         * :rtype: floatType &dpYield
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        druckerPragerSurface(stress, A, B, dpYield);
    
        return dpYield;
    }

    floatType druckerPragerSurface(const floatVector &stress, const floatVector &dpParam){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * :param floatType &vonMises: The von Mises stress
         * :param floatType &meanStress: The mean Stress
         * :param floatVector &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :returns: The Drucker-Prager yield stress/criterion/surface
         * :rtype: floatType &dpYield
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        druckerPragerSurface(stress, dpParam, dpYield);
    
        return dpYield;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :param floatVector &jacobian: The row major jacobian tensor w.r.t. the stress
         */

        //Check vector lengths
        unsigned int length = stress.size();
        if (length != jacobian.size()){
            return new errorNode("druckerPragerSurface", "The stress tensor and jacobian tensor sizes must match.");
        }

        //Calculate von Mises and mean stresses with jacobians
        floatType vonMises, meanStress;
        vonMises = meanStress = 0.;
        floatVector vonMisesJacobian(jacobian.size(), 0.);
        floatVector meanStressJacobian(jacobian.size(), 0.);
        calculateVonMisesStress(stress, vonMises, vonMisesJacobian);
        calculateMeanStress(stress, meanStress, meanStressJacobian);

        //Calculate the Drucker-Prager yield criterion
        druckerPragerSurface(stress, A, B, dpYield); 

        //Calculate the Drucker-Prager jacobian
        jacobian = vonMisesJacobian - A * meanStressJacobian;
    
        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The stress tensor
         * :param floatVector &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :param floatVector &jacobian: The row major jacobian tensor w.r.t. the stress
         */

        //Check Drucker-Prager parameter vector length 
        if (dpParam.size() != 2){
            return new errorNode("druckerPragerSurface", "Two parameters are required for the Drucker-Prager surface.");
        }

        druckerPragerSurface(stress, dpParam[0], dpParam[1], dpYield, jacobian);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :param floatVector &jacobian: The row major jacobian tensor w.r.t. the stress
         * :param floatMatrix &djacobiandstress: The gradient of the jacobian w.r.t. the stress
         */

        //Check vector lengths
        unsigned int length = stress.size();
        if (length != jacobian.size()){
            return new errorNode("druckerPragerSurface", "The stress tensor and jacobian tensor sizes must match.");
        }

        //Calculate von Mises and mean stresses with jacobians
        floatType vonMises, meanStress;
        vonMises = meanStress = 0.;
        floatVector vonMisesJacobian(jacobian.size(), 0.);
        floatVector meanStressJacobian(jacobian.size(), 0.);
        calculateVonMisesStress(stress, vonMises, vonMisesJacobian);
        calculateMeanStress(stress, meanStress, meanStressJacobian);

        //Calculate the Drucker-Prager yield criterion
        druckerPragerSurface(stress, A, B, dpYield); 

        //Calculate the Drucker-Prager jacobian
        jacobian = vonMisesJacobian - A * meanStressJacobian;
    
        //Compute the gradient of the jacobian w.r.t. the stress
        floatVector eye(stress.size(), 0);
        vectorTools::eye<floatType>(eye);
        floatMatrix EYE = vectorTools::eye<floatType>(stress.size());

        floatVector deviatoric = calculateDeviatoricStress(stress);

        djacobiandstress = (3/(2*vonMises))*(EYE - vectorTools::dyadic(eye, meanStressJacobian)
                                                 - vectorTools::dyadic(deviatoric, vonMisesJacobian)/vonMises);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :param floatVector &jacobian: The row major jacobian tensor w.r.t. the stress
         * :param floatVector &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        druckerPragerSurface(stress, A, B, dpYield, jacobian); 

        //Calculate the Drucker-Prager unit normal flow direction as the normalized jacobian
        unitDirection = jacobian / std::sqrt(3./2. + pow(A, 2.)/3.);
    
        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         * f = \sigma^{vonMises} - A*\sigma^{mean} - B
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * :param floatVector &stress: The stress tensor
         * :param floatType &A: The first Drucker-Prager material parameter 
         * :param floatType &B: The second Drucker-Prager material parameter 
         * :param floatType &dpYield: The Drucker-Prager yield stress/criterion/surface
         * :param floatVector &jacobian: The row major jacobian tensor w.r.t. the stress
         * :param floatVector &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * :param floatMatrix &unitDirectionJacobian: The jacobian of the unit direction w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        floatMatrix djacobiandstress;
        druckerPragerSurface(stress, A, B, dpYield, jacobian, djacobiandstress); 

        //Compute the unit normal flow direction and the jacobian of the unit normal flow direction
        //w.r.t. stress
        floatMatrix duDdjacobian;
        constitutiveTools::computeUnitNormal(jacobian, unitDirection, duDdjacobian);

        unitDirectionJacobian = floatMatrix(stress.size(), floatVector(stress.size(), 0));
        for (unsigned int I=0; I<stress.size(); I++){
            for (unsigned int J=0; J<stress.size(); J++){
                for (unsigned int K=0; K<stress.size(); K++){
                    unitDirectionJacobian[I][J] += duDdjacobian[I][K]*djacobiandstress[K][J];
                }
            }
        }
    
        return NULL;
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

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatType n, floatType &p){
        /*!
         * Implementation of the Peryzna type model of the form
         * 
         * p = A \left \langle \frac{f}{q} \right \rangle^n
         * 
         * where \langle \rangle are the Macaulay brackets.
         * 
         * :param const floatType f: The numerator term in the brackets.
         * :param const floatType q: The denominator term in the brackets.
         * :param const floatType A: The scaling factor.
         * :param const floatType n: The exponent.
         * :param const floatType &p: The value of the model.
         */

        if (vectorTools::fuzzyEquals(q, 0.)){
            return new errorNode("peryznaModel", "The denominator term is zero");
        }

        if (n < 1){
            return new errorNode("peryznaModel (jacobian)", "n must be >= 1");
        }

        p = A*pow(constitutiveTools::mac(f/q), n);
        return NULL;
    }

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatType n, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA){

        /*!
         * Implementation of the Peryzna type model of the form
         * 
         * p = A \left \langle \frac{f}{q} \right \rangle^n
         * 
         * where \langle \rangle are the Macaulay brackets.
         * 
         * :param const floatType f: The numerator term in the brackets.
         * :param const floatType q: The denominator term in the brackets.
         * :param const floatType A: The scaling factor.
         * :param const floatType n: The exponent.
         * :param floatType &p: The value of the model.
         * :param floatType &dpdf: The derivative of the value w.r.t. f.
         * :param floatType &dpdq: The derivative of the value w.r.t. q.
         * :param floatType &dpdA: The derivative of the value w.r.t. A.
         */

        if (vectorTools::fuzzyEquals(q, 0.)){
            return new errorNode("peryznaModel (jacobian)", "The denominator term is zero");
        }

        if (n < 1){
            return new errorNode("peryznaModel (jacobian)", "n must be >= 1");
        }

        //Compute the value
        floatType mac, dmacdx;
        mac = constitutiveTools::mac(f/q, dmacdx);

        p = A*pow(constitutiveTools::mac(f/q), n);

        dpdf = A*n*pow(mac, n-1)*dmacdx/q;
        dpdq = -A*n*pow(mac, n-1)*dmacdx*f/(q*q);
        dpdA = pow(constitutiveTools::mac(f/q), n);

        return NULL;
    }
}
