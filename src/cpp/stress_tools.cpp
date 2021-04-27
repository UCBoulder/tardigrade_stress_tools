/**
  ******************************************************************************
  * \file stress_tools.cpp
  ******************************************************************************
  *  A collection of tools which implement and solve stress-strain relationships
  *  in such a way to enable more rapid development of constitutive models which
  *  have capabilities which may not be contained within a collection of
  *  of constitutive models.
  ******************************************************************************
  */

#include<stress_tools.h>

namespace stressTools{

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{1}{3}*trace(\sigma)\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        floatType trace;
        vectorTools::trace(stress, trace);
        meanStress = 1./3.*trace;

        return NULL;
    }

    floatType calculateMeanStress(const floatVector &stress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{1}{3}*trace(\sigma)\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        floatType meanStress;
        calculateMeanStress(stress, meanStress);

        return meanStress;
    }

    errorOut calculateMeanStress(const floatMatrix &stress, floatType &meanStress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         *
         * \f$meanStress = \frac{1}{3}*trace(\sigma)\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        floatType trace;
        vectorTools::trace(stress, trace);
        meanStress = 1./3.*trace;

        return NULL;
    }

    floatType calculateMeanStress(const floatMatrix &stress){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         *
         * \f$meanStress = \frac{1}{3}*trace(\sigma)\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        floatType meanStress;
        calculateMeanStress(stress, meanStress);

        return meanStress;
    }

    errorOut calculateMeanStress(const floatVector &stress, floatType &meanStress, floatVector &jacobian){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{1}{3}*trace(\sigma)\f$
         *
         * \param &stress: The row major stress tensor
         * \param &meanStress: The scalar mean stress
         * \param &jacobian: The row major jacobian w.r.t. the stress tensor
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
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &deviatoric: The deviatoric stress tensor in row major format
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

    errorOut calculateDeviatoricStress(const floatVector &stress, floatVector &deviatoric, floatMatrix &jacobian){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * Also return the jacobian
         *
         * \f$\frac{\partial \sigma_{ij}^{deviatoric}}{\partial \sigma_{kl}} = \delta_{ik}\delta{jl} - \frac{\partial
         * \bar{\sigma}}{\partial \sigma_{kl}} \delta_{ij}\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &deviatoric: The deviatoric stress tensor in row major format
         * \param &jacobian: The jacobian of the deviatoric stress tensor w.r.t. the stress.
         */

        //Compute the deviatoric stress.
        deviatoric.resize(stress.size());
        errorOut error = calculateDeviatoricStress(stress, deviatoric);

        if (error){
            errorOut result = new errorNode("calculateDeviatoricStress (jacobian) ", "Error in calculation of deviatoric stress");
            result->addNext(error);
            return result;
        }

        //Compute the jacobian
        floatVector eye(stress.size(), 0);
        vectorTools::eye(eye);

        jacobian = vectorTools::eye<floatType>(stress.size()) - 1./3 * vectorTools::dyadic(eye, eye);

        return NULL;
    }

    floatVector calculateDeviatoricStress(const floatVector &stress){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * \param &stress: The stress tensor in row major format
         * \returns <B>deviatoric<\B>: The deviatoric stress tensor in row major format
         */

        floatVector deviatoric(stress.size());
        errorOut result = calculateDeviatoricStress(stress, deviatoric);
        if (result){
            result->print();
            throw std::runtime_error("Error in calculation of deviatoric stress");
        }

        return deviatoric;
    }

    floatVector calculateDeviatoricStress(const floatVector &stress, floatMatrix &jacobian){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * Also return the jacobian
         *
         * \f$\frac{\partial \sigma_{ij}^{deviatoric}}{\partial \sigma_{kl}} = \delta_{ik}\delta{jl} - \frac{\partial
         * \bar{\sigma}}{\partial \sigma_{kl}} \delta_{ij}\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &jacobian: The jacobian of the deviatoric stress tensor w.r.t. the stress.
         */

        floatVector deviatoric(stress.size());
        errorOut result = calculateDeviatoricStress(stress, deviatoric, jacobian);
        if (result){
            result->print();
            throw std::runtime_error("Error in calculation of deviatoric stress");
        }

        return deviatoric;
    }

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}\f$
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * \param &stress: The row major stress tensor
         * \param &meanStress: The scalar mean stress
         */

        floatVector deviatoric = calculateDeviatoricStress(stress);
        vonMises = std::sqrt(3./2.*vectorTools::inner(deviatoric, deviatoric));

        return NULL;
    }

    floatType calculateVonMisesStress(const floatVector &stress){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}\f$
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * \param &stress: The row major stress tensor
         * \param &meanStress: The scalar mean stress
         */

        floatType vonMises = 0.;
        errorOut result = calculateVonMisesStress(stress, vonMises);
        if (result){
            result->print();
            throw std::runtime_error("Error in calculation of deviatoric stress");
        }

        return vonMises;
    }

    errorOut calculateVonMisesStress(const floatVector &stress, floatType &vonMises, floatVector &jacobian){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{vonMises} = \sqrt{\frac{3}{2}*\sigma^{deviatoric}\sigma^{deviatoric}}\f$
         *
         * \f$\sigma^{deviatoric} = \sigma - \sigma^{mean}I\f$
         *
         * \param &stress: The row major stress tensor
         * \param &meanStress: The scalar mean stress
         * \param &jacobian: The row major mean stress jacobian tensor w.r.t. the stress tensor
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
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        dpYield = vonMises + A*meanStress - B;

        return NULL;
    }

    errorOut druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
         *
         * \f$f = \sigma^{vonMises} + dpParam[0]*\sigma^{mean} - dpParam[1]\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
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
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \returns <B>dpYield<\B>: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0;

        //Calculate DP yield criterion
        druckerPragerSurface(vonMises, meanStress, A, B, dpYield);

        return dpYield;
    }

    floatType druckerPragerSurface(const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
         *
         * \f$f = \sigma^{vonMises} + dpParam[0]*\sigma^{mean} - dpParam[1]\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \returns <B>dpYield<\B>: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0;
        druckerPragerSurface(vonMises, meanStress, dpParam, dpYield);

        return dpYield;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The row major stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
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
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
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
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The row major stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \returns <B>dpYield</B>: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        druckerPragerSurface(stress, A, B, dpYield);

        return dpYield;
    }

    floatType druckerPragerSurface(const floatVector &stress, const floatVector &dpParam){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \returns <B>dpYield<\B>: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        druckerPragerSurface(stress, dpParam, dpYield);

        return dpYield;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
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
        jacobian = vonMisesJacobian + A * meanStressJacobian;

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
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
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &djacobiandstress: The gradient of the jacobian w.r.t. the stress
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
        jacobian = vonMisesJacobian + A * meanStressJacobian;

        //Compute the gradient of the jacobian w.r.t. the stress
        floatVector eye(stress.size(), 0);
        vectorTools::eye<floatType>(eye);
        floatMatrix EYE = vectorTools::eye<floatType>(stress.size());

        floatVector deviatoric = calculateDeviatoricStress(stress);

        djacobiandstress = (3/(2*vonMises))*(EYE - vectorTools::dyadic(eye, meanStressJacobian)
                                                 - vectorTools::dyadic(deviatoric, vonMisesJacobian)/vonMises);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &djacobiandstress: The gradient of the jacobian w.r.t. the stress
         */

        //Check Drucker-Prager parameter vector length
        if (dpParam.size() != 2){
            return new errorNode("druckerPragerSurface", "Two parameters are required for the Drucker-Prager surface.");
        }

        druckerPragerSurface(stress, dpParam[0], dpParam[1], dpYield, jacobian, djacobiandstress);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        druckerPragerSurface(stress, A, B, dpYield, jacobian);

        //Calculate the Drucker-Prager unit normal flow direction as the normalized jacobian
        unitDirection = jacobian / std::sqrt(3./2. + pow(A, 2.)/3.);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        druckerPragerSurface(stress, dpParam, dpYield, jacobian);

        //Calculate the Drucker-Prager unit normal flow direction as the normalized jacobian
        unitDirection = jacobian / std::sqrt(3./2. + pow(dpParam[0], 2.)/3.);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * \param &unitDirectionJacobian: The jacobian of the unit direction w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        floatMatrix djacobiandstress;
        druckerPragerSurface(stress, A, B, dpYield, jacobian, djacobiandstress);

        //Compute the unit normal flow direction and the jacobian of the unit normal flow direction
        //w.r.t. stress
        floatMatrix duDdjacobian;
        constitutiveTools::computeUnitNormal(jacobian, unitDirection, duDdjacobian);

        unitDirectionJacobian = vectorTools::dot(duDdjacobian, djacobiandstress);

        return NULL;
    }

    errorOut druckerPragerSurface(const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{vonMises} + A*\sigma^{mean} - B\f$
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector {A, B}
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * \param &unitDirectionJacobian: The jacobian of the unit direction w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        floatMatrix djacobiandstress;
        druckerPragerSurface(stress, dpParam, dpYield, jacobian, djacobiandstress);

        //Compute the unit normal flow direction and the jacobian of the unit normal flow direction
        //w.r.t. stress
        floatMatrix duDdjacobian;
        constitutiveTools::computeUnitNormal(jacobian, unitDirection, duDdjacobian);

        unitDirectionJacobian = vectorTools::dot(duDdjacobian, djacobiandstress);

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
         *
         * \f$\rho^0 \Psi = 0.5*(E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ}))\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * \param &previousStrain: The previous value of strain
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [\f$G_{\infty}\f$, \f$\tau\f$ s, \f$G\f$ s] where
         *         - \f$G_{\infty}\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter (0 for implicit, 1 for explicit)
         * \param &stress: The computed stress in the reference configuration (i.e. the same configuration as the strain)
         * \param &currentStateVariables: The current values of the state variables
         */

        floatType dt = currentTime - previousTime;

        //Check the material parameters
        if (materialParameters.size() == 0){
            return new errorNode("linearViscoelasticity", "At least one material parameter must be provided.");
        }

        //Set the number of Prony-Series terms
        if ( (materialParameters.size() - 1) % 2 != 0 ){
            return new errorNode("linearViscoelasticity", "An equal number of taus and Gs must be provided.");
        }
        unsigned int nTerms = (materialParameters.size() - 1)/2;

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
         *
         * \f$\rho^0 \Psi = 0.5*(E_{IJ} G_{\infty} E_{IJ} + \sum_{n=1}^N (E_{IJ} - \Xi_{IJ}^n) G^n (E_{IJ} - \Xi_{IJ}))\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &previousStrain: The previous value of strain
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or (potentially) other non-linear effects.
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [\f$G_{\infty}\f$, \f$\tau\f$ s, \f$G\f$ s] where
         *         - \f$G_{\infty}\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter (0 for implicit, 1 for explicit)
         * \param &stress: The computed stress in the reference configuration (i.e. the same configuration as the strain)
         * \param &currentStateVariables: The current values of the state variables
         * \param &dstressdstrain: The derivative of the stress w.r.t. the strain
         * \param &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
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
         *
         * \f$U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &jacobian: The jacobian of deformation
         * \param &bulkModulus: The bulk modulus
         * \param NOTE: It is up to the user to determine
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
         *
         * \f$U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &jacobian: The jacobian of deformation
         * \param &bulkModulus: The bulk modulus
         * \param NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         * \param &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian
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
         *
         * \f$U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &bulkModulus: The bulk modulus
         * \param NOTE: It is up to the user to determine
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
         *
         * \f$U(J) = 0.5*bulkModulus*(0.5*(J**2 - 1) - ln(J))\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &bulkModulus: The bulk modulus
         * \param NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         * \param &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian
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
         * \f$p = A \left \langle \frac{f}{q} \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param n: The exponent.
         * \param &p: The value of the model.
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

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p){
        /*!
         * Implementation of the Peryzna type model of the form
         *
         * \f$p = A \left \langle \frac{f}{q} \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param &parameters: The parameters (currently just n)
         * \param &p: The value of the model.
         */
        if (parameters.size() != 1){
            return new errorNode("peryznaModel", "The parameters vector is one value long");
        }
        return peryznaModel(f, q, A, parameters[0], p);
    }

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatType n, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA){

        /*!
         * Implementation of the Peryzna type model of the form
         *
         * \f$p = A \left \langle \frac{f}{q} \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param n: The exponent.
         * \param &p: The value of the model.
         * \param &dpdf: The derivative of the value w.r.t. f.
         * \param &dpdq: The derivative of the value w.r.t. q.
         * \param &dpdA: The derivative of the value w.r.t. A.
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

    errorOut peryznaModel(const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA){
        /*!
         * Implementation of the Peryzna type model of the form
         *
         * \f$p = A \left \langle \frac{f}{q} \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param &parameters: The parameters (currently just n)
         * \param &p: The value of the model.
         * \param &dpdf: The derivative of the value w.r.t. f.
         * \param &dpdq: The derivative of the value w.r.t. q.
         * \param &dpdA: The derivative of the value w.r.t. A.
         */
        if (parameters.size() != 1){
            return new errorNode("peryznaModel", "The parameters vector is one value long");
        }
        return peryznaModel(f, q, A, parameters[0], p, dpdf, dpdq, dpdA);
    }

    errorOut linearHardening(const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value){
        /*!
         * Compute the linear hardening curve value.
         *
         * \f$value = stateVariables_i linearModuli_i + scalarShift\f$
         *
         * \param &stateVariables: The state variable vector
         * \param &linearModuli: The linear moduli vector.
         * \param &scalarShift: The scalar shift value.
         * \param &value: The value of the linear hardening curve.
         */

        if (stateVariables.size() != linearModuli.size()){
            return new errorNode("linearHardening", "The state variables and the moduli must have the same size");
        }

        value = vectorTools::dot(stateVariables, linearModuli) + scalarShift;
        return NULL;
    }

    errorOut linearHardening(const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value, floatVector &valueJacobian){
        /*!
         * Compute the linear hardening curve value.
         *
         * \f$value = stateVariables_i linearModuli_i + scalarShift\f$
         *
         * \param &stateVariables: The state variable vector
         * \param &linearModuli: The linear moduli vector.
         * \param &scalarShift: The scalar shift value.
         * \param &value: The value of the linear hardening curve.
         * \param &valueJacobian: The jacobian of value w.r.t. the state variables.
         */

        errorOut error = linearHardening(stateVariables, linearModuli, scalarShift, value);

        if (error){
            errorOut result = new errorNode("linearHardening (jacobian)", "Error in computation of hardening value");
            result->addNext(error);
            return result;
        }

        valueJacobian = linearModuli;
        return NULL;
    }

    errorOut computeJaumannStiffnessTensor( const floatVector &cauchyStress, const floatVector &currentDeformationGradient,
                                            const floatMatrix &dCauchydF, floatMatrix &C ){
        /*!
         * Compute the Jaumann stiffness tensor from the cauchy stress, the current deformation gradient,
         * and the total derivative of the Cauchy stress w.r.t. the deformation gradient
         * 
         * From the variation of the Jaumann rate
         * \f$J \mathbb{C}_{ijkl} \delta D_{kl} = \delta \left(J \sigma_{ij}\right) + \delta W_{ik} \sigma_{kj} - \sigma_{ik} \delta W_{kj}\f$
         * 
         * Where \f$J$\f is the Jacobian of deformation, \f$\mathbb{C}\f$ is the Jaumann stiffness tensor,
         * \f$\bf{\sigma}\f$ is the Cauchy stress and \f$\delta \bf{D}\f$ is the perturbation of the rate of
         * deformation, and \f$\delta \bf{W}\f$ is the perturbation of the rate of spin.
         * 
         * By using the properties that
         * \f$\delta D_{kl} = \text{symm}\left(\delta F_{kK} F_{Kl}^{-1}\right)\f$
         * \f$\delta W_{kl} = \text{asymm}\left(\delta F_{kK} F_{Kl}^{-1}\right)\f$
         * 
         * Where
         * \f$\text{symm}\left(\bf{A}\right) = \frac{1}{2}\left(\bf{A} + \bf{A}^T\right)\f$
         * \f$\text{asymm}\left(\bf{A}\right) = \frac{1}{2}\left(\bf{A} - \bf{A}^T\right)\f$
         * 
         * It can be shown that
         * \f$\mathbb{C}_{ijkl} = \delta_{kl} \sigma_{ij} + \frac{D \sigma_{ij}}{D F_{kK}} F_{lK} + \mathbb{P}_{irkl}^{asymm} \sigma_{rj} - \mathbb{P}_{rjkl}^{asymm} \sigma_{ir}\f$
         * 
         * Where
         * \f$\mathbb{P}_{ijkl}^{asymm} = \frac{1}{2}\left(\delta_{ik} \delta_{jl} - \delta_{jk} \delta_{il}\right)\f$
         * 
         * \param &cauchyStress: The Cauchy stress
         * \param &currentDeformationGradient: The current deformation gradient i.e. the mapping for differential
         *     lengths from the reference to the current configuration
         * \param &dCauchydF: The Jacobian (total derivative) of the Cauchy stress w.r.t. the deformation
         *     gradient
         * \param &C: The Jaumann stiffness tensor
         */

        // Build the second order identity tensor
        floatVector eye( cauchyStress.size( ), 0 );
        vectorTools::eye( eye );

        // Set the dimension of the problem
        unsigned int dim = vectorTools::trace( eye );

        // Construct the asymmetric projection tensor as a vector
        floatVector Pasymm( cauchyStress.size( ) * cauchyStress.size( ), 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    for ( unsigned int l = 0; l < dim; l++ ){

                        Pasymm[ dim * dim * dim * i + dim * dim * j + dim * k + l ] =
                            eye[ dim * i + k ] * eye[ dim * j + l ] - eye[ dim * j + k ] * eye[ dim * i + l ];

                    }

                }

            }

        }

        // Initialize the stiffness tensor by including cauchyStress dyad eye
        C = vectorTools::dyadic( cauchyStress, eye );

        // Add the dCauchydF and skew symmetric terms

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    for ( unsigned int l = 0; l < dim; l++ ){

                        for ( unsigned int r = 0; r < dim; r++ ){

                            C[ dim * dim * dim * i + dim * dim * j + dim * k + l ] +=
                                dCauchydF[ dim * i + j ][ dim * k + r ] * currentDeformationGradient[ dim * l + r ]
                              + Pasymm[ dim * dim * dim * i + dim * dim * r + dim * k + l ] * cauchyStress[ dim * r + j ]
                              - Pasymm[ dim * dim * dim * r + dim * dim * j + dim * k + l ] * cauchyStress[ dim * i + r ];

                        }

                    }

                }

            }

        }

        return NULL;

    }

}
