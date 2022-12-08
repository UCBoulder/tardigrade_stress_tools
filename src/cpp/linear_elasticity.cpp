#include<linear_elasticity.h>
#include<sstream>

namespace linearElasticity{

    /** Define the expected number of spatial dimensions */
    unsigned int spatialDimensions = 3;

    errorOut formReferenceStiffnessTensor( const floatVector &parameters, floatMatrix &C ){
        /*!
         * Form the stiffness tensor in the reference configuration.
         * 
         * \param &parameters: The parameters. The first index is the type of stiffness tensor and the later values are the coefficients.
         *     -type 0: Isotropic stiffness parameterized by lambda and mu
         * \param &C: The resulting stiffness tensor.
         */

        if ( parameters.size( ) < 1 ){

            return new errorNode( __func__, "The parameters must at least define the type" );

        }

        // Form isotropic stiffness tensor
        if ( int( parameters[ 0 ] + 0.5 ) == 0 ){

            if ( parameters.size( ) != 3 ){

                return new errorNode( __func__, "Isotropic stiffeness requires two parameters. Parameters only defines " + std::to_string( parameters.size( ) - 1 ) );

            }

            floatType lambda = parameters[ 1 ];
            floatType mu     = parameters[ 2 ];

            C = {
                    { lambda + 2 * mu,      0,      0,      0,          lambda,      0,      0,      0,          lambda },
                    {               0,      0,      0, 2 * mu,               0,      0,      0,      0,               0 },
                    {               0,      0,      0,      0,               0,      0, 2 * mu,      0,               0 },
                    {               0, 2 * mu,      0,      0,               0,      0,      0,      0,               0 },
                    {          lambda,      0,      0,      0, lambda + 2 * mu,      0,      0,      0,          lambda },
                    {               0,      0,      0,      0,               0,      0,      0, 2 * mu,               0 },
                    {               0,      0, 2 * mu,      0,               0,      0,      0,      0,               0 },
                    {               0,      0,      0,      0,               0, 2 * mu,      0,      0,               0 },
                    {          lambda,      0,      0,      0,          lambda,      0,      0,      0, lambda + 2 * mu }
                }; 

        }
        else{

            return new errorNode( __func__, "Type " + std::to_string( parameters[ 0 ] ) + " is not recognized" );

        }

        return NULL;

    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] \f$
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {
            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        floatMatrix C;

        error = formReferenceStiffnessTensor( parameters, C );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the reference stiffness tensor" );
            result->addNext( error );
            return result;

        }

        energy = 0.5 * vectorTools::dot( vectorTools::dot( C, E ), E ) / detChi;

        return NULL;
    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] \f$
         * 
         * and the value of the Cauchy stress which is defined via
         * 
         * \f$\sigma_{ij} = \frac{1}{J} \frac{ \partial \left( \rho \psi \right ) }{\partial F_{iI}} F_{jI} \f$
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         * \param &cauchyStress; The expected cauchy stress
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        floatMatrix dEdChi;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E, dEdChi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;
        floatVector dDetChidChi;
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {

            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        floatMatrix C;
        error = formReferenceStiffnessTensor( parameters, C );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the reference stiffness tensor" );
            result->addNext( error );
            return result;

        }

        floatVector eye( spatialDimensions * spatialDimensions );
        vectorTools::eye( eye );

        energy = 0.5 * vectorTools::dot( vectorTools::dot( C, E ), E ) / detChi;

        floatVector dEnergydChi = vectorTools::Tdot( dEdChi, vectorTools::dot( C, E ) ) / detChi
                                - energy * vectorTools::matrixMultiply( eye, invChi, spatialDimensions, spatialDimensions, spatialDimensions, spatialDimensions, false, true );

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;
                }
            }
        }


        return NULL;
    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] \f$
         * 
         * and the value of the Cauchy stress which is defined via
         * 
         * \f$\sigma_{ij} = \frac{1}{J} \frac{ \partial \left( \rho \psi \right ) }{\partial F_{iI}} F_{jI} \f$
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         * \param &cauchyStress; The expected cauchy stress
         * \param &dEnergydChi: The gradient of the energy w.r.t. the micro deformation
         * \param &dCauchyStressdChi: The gradient of the Cauchy stress w.r.t. the micro deformation
         * 
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        floatMatrix dEdChi;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E, dEdChi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;
        floatVector dDetChidChi;
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {

            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        floatMatrix C;
        error = formReferenceStiffnessTensor( parameters, C );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the reference stiffness tensor" );
            result->addNext( error );
            return result;

        }

        floatVector eye( spatialDimensions * spatialDimensions );
        vectorTools::eye( eye );

        floatVector invChiT = vectorTools::matrixMultiply( eye, invChi, spatialDimensions, spatialDimensions, spatialDimensions, spatialDimensions, false, true );

        floatVector CE = vectorTools::dot( C, E );
        floatMatrix dCEdChi = vectorTools::dot( C, dEdChi );

        energy = 0.5 * vectorTools::dot( CE, E ) / detChi;

        dEnergydChi = vectorTools::Tdot( dEdChi, CE ) / detChi - energy * invChiT;

        floatVector d2EnergydChi2( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0 );

        for ( unsigned int k = 0; k < spatialDimensions; k++ ){
            for ( unsigned int K = 0; K < spatialDimensions; K++ ){
                for ( unsigned int l = 0; l < spatialDimensions; l++ ){
                    for ( unsigned int L = 0; L < spatialDimensions; L++ ){
                        d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                            += eye[ spatialDimensions * l + k ] * CE[ spatialDimensions * K + L ] / detChi
                             - dEnergydChi[ spatialDimensions * l + L ] * invChi[ spatialDimensions * K + k ]
                             + energy * invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + k ];
                        for ( unsigned int IJ = 0; IJ < spatialDimensions * spatialDimensions; IJ++ ){
                            d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                                += ( dEdChi[ IJ ][ spatialDimensions * k + K ] * dCEdChi[ IJ ][ spatialDimensions * l + L ]
                                   - dEdChi[ IJ ][ spatialDimensions * k + K ] * CE[ IJ ] * invChi[ spatialDimensions * L + l ] ) / detChi;
                        }
                    }
                }
            }
        }

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        dCauchyStressdChi = floatMatrix( spatialDimensions * spatialDimensions, floatVector( spatialDimensions * spatialDimensions, 0 ) );

        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;

                    for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                        for ( unsigned int K = 0; K < spatialDimensions; K++ ){

                            dCauchyStressdChi[ spatialDimensions * i + j ][ spatialDimensions * k + K ] += ( d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * k + K ] * chi[ spatialDimensions * j + I ]
                                                                                                           + dEnergydChi[ spatialDimensions * i + I ] * eye[ spatialDimensions * j + k ] * eye[ spatialDimensions * I + K ]
                                                                                                           - dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] * invChi[ spatialDimensions * K + k ] ) / detChi;

                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut evaluateEnergy( const floatVector &chi, const floatVector &parameters, floatType &energy, floatVector &cauchyStress,
                             floatVector &dEnergydChi, floatMatrix &dCauchyStressdChi,
                             floatVector &d2EnergydChi2, floatMatrix &d2CauchyStressdChi2 ){
        /*!
         * Compute the value of the linear elastic energy which we define via
         * 
         * \f$\rho \psi = \frac{1}{J} \left[ \frac{\lambda}{2} \left( E_{II} \right)^2 + \mu E_{IJ} E_{JI} \right] \f$
         * 
         * and the value of the Cauchy stress which is defined via
         * 
         * \f$\sigma_{ij} = \frac{1}{J} \frac{ \partial \left( \rho \psi \right ) }{\partial F_{iI}} F_{jI} \f$
         * 
         * \param &chi: The micro-deformation
         * \param &parameters: The parameters used in the calculation. The two Lame parameters are expected lambda and mu.
         * \param &energy: The resulting free energy in the current configuration
         * \param &cauchyStress; The expected cauchy stress
         * \param &dEnergydChi: The gradient of the energy w.r.t. the micro deformation
         * \param &dCauchyStressdChi: The gradient of the Cauchy stress w.r.t. the micro deformation
         * \param &d2EnergydChi2: The second gradient of the energy w.r.t. the micro deformation
         * \param &d2CauchyStressdChi2: The second gradient of the Cauchy stress w.r.t. the micro deformation
         * 
         */

        if ( chi.size( ) != spatialDimensions * spatialDimensions ){
            return new errorNode( __func__, "The spatial dimension of " + std::to_string( spatialDimensions ) + " is not reflected by chi which has a size of " + std::to_string( chi.size( ) ) + " rather than " + std::to_string( spatialDimensions * spatialDimensions ) );
        }

        floatVector E;
        floatMatrix dEdChi;
        errorOut error = constitutiveTools::computeGreenLagrangeStrain( chi, E, dEdChi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        floatType detChi;
        floatVector dDetChidChi;
        floatVector invChi;

        try{

            detChi = vectorTools::determinant( chi, spatialDimensions, spatialDimensions );
            dDetChidChi = vectorTools::computeDDetAdJ( chi, spatialDimensions, spatialDimensions );

            invChi = vectorTools::inverse( chi, spatialDimensions, spatialDimensions );

        }
        catch ( std::exception &e ) {

            std::ostringstream message;
            message << "Error in calculation of det( chi ). Error follows:\n";
            message << e.what( );
            return new errorNode( __func__, message.str( ) );

        }

        floatMatrix C;
        error = formReferenceStiffnessTensor( parameters, C );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the reference stiffness tensor" );
            result->addNext( error );
            return result;

        }

        floatVector eye( spatialDimensions * spatialDimensions );
        vectorTools::eye( eye );

        floatVector invChiT = vectorTools::matrixMultiply( eye, invChi, spatialDimensions, spatialDimensions, spatialDimensions, spatialDimensions, false, true );

        floatVector CE = vectorTools::dot( C, E );
        floatMatrix dCEdChi = vectorTools::dot( C, dEdChi );

        energy = 0.5 * vectorTools::dot( CE, E ) / detChi;

        dEnergydChi = vectorTools::Tdot( dEdChi, CE ) / detChi - energy * invChiT;

        d2EnergydChi2 = floatVector( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0 );

        floatVector d3EnergydChi3( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0 );

        for ( unsigned int k = 0; k < spatialDimensions; k++ ){
            for ( unsigned int K = 0; K < spatialDimensions; K++ ){
                for ( unsigned int l = 0; l < spatialDimensions; l++ ){
                    for ( unsigned int L = 0; L < spatialDimensions; L++ ){
                        d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                            += eye[ spatialDimensions * l + k ] * CE[ spatialDimensions * K + L ] / detChi
                             - dEnergydChi[ spatialDimensions * l + L ] * invChi[ spatialDimensions * K + k ]
                             + energy * invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + k ];
                        for ( unsigned int IJ = 0; IJ < spatialDimensions * spatialDimensions; IJ++ ){
                            d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                                += ( dEdChi[ IJ ][ spatialDimensions * k + K ] * dCEdChi[ IJ ][ spatialDimensions * l + L ]
                                   - dEdChi[ IJ ][ spatialDimensions * k + K ] * CE[ IJ ] * invChi[ spatialDimensions * L + l ] ) / detChi;
                        }
                    }
                }
            }
        }
        for ( unsigned int k = 0; k < spatialDimensions; k++ ){
            for ( unsigned int K = 0; K < spatialDimensions; K++ ){
                for ( unsigned int l = 0; l < spatialDimensions; l++ ){
                    for ( unsigned int L = 0; L < spatialDimensions; L++ ){
                        for ( unsigned int m = 0; m < spatialDimensions; m++ ){
                            for ( unsigned int M = 0; M < spatialDimensions; M++ ){
                                d3EnergydChi3[ spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * K + spatialDimensions * spatialDimensions * spatialDimensions * l + spatialDimensions * spatialDimensions * L + spatialDimensions * m + M ]
                                    += ( eye[ spatialDimensions * l + k ] * dCEdChi[ spatialDimensions * K + L ][ spatialDimensions * m + M ] - eye[ spatialDimensions * l + k ] * CE[ spatialDimensions * K + L ] * invChi[ spatialDimensions * M + m ] ) / detChi
                                     - d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * l + spatialDimensions * spatialDimensions * L + spatialDimensions * m + M ] * invChi[ spatialDimensions * K + k ]
                                     + dEnergydChi[ spatialDimensions * l + L ] * invChi[ spatialDimensions * K + m ] * invChi[ spatialDimensions * M + k ]
                                     + dEnergydChi[ spatialDimensions * m + M ] * invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + k ]
                                     - energy * invChi[ spatialDimensions * K + m ] * invChi[ spatialDimensions * M + l ] * invChi[ spatialDimensions * L + k ]
                                     - energy * invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + m ] * invChi[ spatialDimensions * M + k ]
                                     + ( eye[ spatialDimensions * m + k ] * dCEdChi[ spatialDimensions * K + M ][ spatialDimensions * l + L ] 
                                       - eye[ spatialDimensions * m + k ] * CE[ spatialDimensions * K + M ] * invChi[ spatialDimensions * L + l ]
                                       )/ detChi;
                                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                                    for ( unsigned int J = 0; J < spatialDimensions; J++ ){
                                        d3EnergydChi3[ spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * K + spatialDimensions * spatialDimensions * spatialDimensions * l + spatialDimensions * spatialDimensions * L + spatialDimensions * m + M ]
                                            += (
                                               + dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * C[ spatialDimensions * I + J ][ spatialDimensions * L + M ] * eye[ spatialDimensions * m + l ]
                                               - dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * dCEdChi[ spatialDimensions * I + J ][ spatialDimensions * m + M ] * invChi[ spatialDimensions * L + l ]
                                               + dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * CE[ spatialDimensions * I + J ] * invChi[ spatialDimensions * L + m ] * invChi[ spatialDimensions * M + l ]
                                               - dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * dCEdChi[ spatialDimensions * I + J ][ spatialDimensions * l + L ] * invChi[ spatialDimensions * M + m ]
                                               + dEdChi[ spatialDimensions * I + J ][ spatialDimensions * k + K ] * CE[ spatialDimensions * I + J ] * invChi[ spatialDimensions * L + l ] * invChi[ spatialDimensions * M + m ]
                                               ) / detChi;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Use the energy gradient to compute the Cauchy stress
        cauchyStress = floatVector( spatialDimensions * spatialDimensions, 0 );
        dCauchyStressdChi = floatMatrix( spatialDimensions * spatialDimensions, floatVector( spatialDimensions * spatialDimensions, 0 ) );
        d2CauchyStressdChi2 = floatMatrix( spatialDimensions * spatialDimensions, floatVector( spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions, 0 ) );

        for ( unsigned int i = 0; i < spatialDimensions; i++ ){
            for ( unsigned int j = 0; j < spatialDimensions; j++ ){
                for ( unsigned int I = 0; I < spatialDimensions; I++ ){
                    cauchyStress[ spatialDimensions * i + j ] += dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] / detChi;

                    for ( unsigned int k = 0; k < spatialDimensions; k++ ){

                        for ( unsigned int K = 0; K < spatialDimensions; K++ ){


                            dCauchyStressdChi[ spatialDimensions * i + j ][ spatialDimensions * k + K ] += ( d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * k + K ] * chi[ spatialDimensions * j + I ]
                                                                                                           + dEnergydChi[ spatialDimensions * i + I ] * eye[ spatialDimensions * j + k ] * eye[ spatialDimensions * I + K ]
                                                                                                           - dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] * invChi[ spatialDimensions * K + k ] ) / detChi;
                            for ( unsigned int l = 0; l < spatialDimensions; l++ ){
                                
                                for ( unsigned int L = 0; L < spatialDimensions; L++ ){

                                    d2CauchyStressdChi2[ spatialDimensions * i + j ][ spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ]
                                        += ( d3EnergydChi3[ spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * spatialDimensions * spatialDimensions * I + spatialDimensions * spatialDimensions * spatialDimensions * k + spatialDimensions * spatialDimensions * K + spatialDimensions * l + L ] * chi[ spatialDimensions * j + I ]
                                           + d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * k + K ] * eye[ spatialDimensions * j + l ] * eye[ spatialDimensions * I + L ]
                                           + d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * l + L ] * eye[ spatialDimensions * j + k ] * eye[ spatialDimensions * I + K ]
                                           - d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * l + L ] * chi[ spatialDimensions * j + I ] * invChi[ spatialDimensions * K + k ]
                                           - dEnergydChi[ spatialDimensions * i + I ] * eye[ spatialDimensions * j + l ] * eye[ spatialDimensions * I + L ] * invChi[ spatialDimensions * K + k ]
                                           + dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] * invChi[ spatialDimensions * K + l ] * invChi[ spatialDimensions * L + k ]
                                           - ( d2EnergydChi2[ spatialDimensions * spatialDimensions * spatialDimensions * i + spatialDimensions * spatialDimensions * I + spatialDimensions * k + K ] * chi[ spatialDimensions * j + I ]
                                             + dEnergydChi[ spatialDimensions * i + I ] * eye[ spatialDimensions * j + k ] * eye[ spatialDimensions * I + K ]
                                             - dEnergydChi[ spatialDimensions * i + I ] * chi[ spatialDimensions * j + I ] * invChi[ spatialDimensions * K + k ]
                                             ) * invChi[ spatialDimensions * L + l ]
                                           ) / detChi;

                                }

                            }

                        }

                    }

                }

            }

        }

        return NULL;

    }

}

