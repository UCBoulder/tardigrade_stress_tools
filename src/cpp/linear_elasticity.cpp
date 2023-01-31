#include<linear_elasticity.h>
#include<sstream>

namespace stressTools{
namespace linearElasticity{

    /** Define the expected number of spatial dimensions */
    unsigned int spatialDimensions = 3;

    errorOut formReferenceStiffnessTensor( const floatVector &parameters, floatMatrix &stiffnessTensor ){
        /*!
         * Form the stiffness tensor in the reference configuration.
         *
         * \f$C_{ijkl} = \begin{bmatrix}
         *   C_{1111} & C_{1112} & C_{1113} & C_{1112} & C_{1122} & C_{1123} & C_{1113} & C_{1123} & C_{1133} \\
         *   C_{1112} & C_{1212} & C_{1213} & C_{1212} & C_{1222} & C_{1223} & C_{1213} & C_{1223} & C_{1233} \\
         *   C_{1113} & C_{1213} & C_{1313} & C_{1213} & C_{1322} & C_{1323} & C_{1313} & C_{1323} & C_{1333} \\
         *   C_{1112} & C_{1212} & C_{1213} & C_{1212} & C_{1222} & C_{1223} & C_{1213} & C_{1223} & C_{1233} \\
         *   C_{1122} & C_{1222} & C_{1322} & C_{1222} & C_{2222} & C_{2223} & C_{1322} & C_{2223} & C_{2233} \\
         *   C_{1123} & C_{1223} & C_{1323} & C_{1223} & C_{2223} & C_{2323} & C_{1323} & C_{2323} & C_{2333} \\
         *   C_{1113} & C_{1213} & C_{1313} & C_{1213} & C_{1322} & C_{1323} & C_{1313} & C_{1323} & C_{1333} \\
         *   C_{1123} & C_{1223} & C_{1323} & C_{1223} & C_{2223} & C_{2323} & C_{1323} & C_{2323} & C_{2333} \\
         *   C_{1133} & C_{1233} & C_{1333} & C_{1233} & C_{2233} & C_{2333} & C_{1333} & C_{2333} & C_{3333}
         * \end{bmatrix}\f$
         *
         * \param &parameters: The tensor components of the 9x9 stiffness tensor. Vector length determines the symmetry.
         *
         * - 81: A row-major vector representing the full 9x9 stiffness tensor directly.
         * - 21: fully anistropic \f$C_{1111}\f$, \f$C_{1112}\f$, \f$C_{1113}\f$, \f$C_{1122}\f$, \f$C_{1123}\f$,
         *   \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1213}\f$, \f$C_{1222}\f$, \f$C_{1223}\f$, \f$C_{1233}\f$,
         *   \f$C_{1313}\f$, \f$C_{1322}\f$, \f$C_{1323}\f$, \f$C_{1333}\f$, \f$C_{2222}\f$, \f$C_{2223}\f$,
         *   \f$C_{2233}\f$, \f$C_{2323}\f$, \f$C_{2333}\f$, \f$C_{3333}\f$
         * - 9: orthotropic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1313}\f$,
         *   \f$C_{2222}\f$, \f$C_{2323}\f$, \f$C_{3333}\f$
         * - 5: transversly isotropic or hexagonal \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1313}\f$, \f$C_{3333}\f$
         * - 3: cubic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1212}\f$
         * - 2: isotropic: lambda (\f$C_{1122}\f$), mu (\f$C_{1212}\f$)
         *
         * \param &stiffnessTensor: The resulting stiffness tensor.
         */
        floatType C1111 = 0.;
        floatType C1112 = 0.;
        floatType C1113 = 0.;
        floatType C1122 = 0.;
        floatType C1123 = 0.;
        floatType C1133 = 0.;
        floatType C1212 = 0.;
        floatType C1213 = 0.;
        floatType C1222 = 0.;
        floatType C1223 = 0.;
        floatType C1233 = 0.;
        floatType C1313 = 0.;
        floatType C1322 = 0.;
        floatType C1323 = 0.;
        floatType C1333 = 0.;
        floatType C2222 = 0.;
        floatType C2223 = 0.;
        floatType C2233 = 0.;
        floatType C2323 = 0.;
        floatType C2333 = 0.;
        floatType C3333 = 0.;

        if ( parameters.size( ) == 81 ){

            unsigned int length = 9;
            stiffnessTensor = vectorTools::inflate( parameters, length, length );
            return NULL;

        }
        else if ( parameters.size( ) == 21 ){

            C1111 = parameters[  0 ];
            C1112 = parameters[  1 ];
            C1113 = parameters[  2 ];
            C1122 = parameters[  3 ];
            C1123 = parameters[  4 ];
            C1133 = parameters[  5 ];
            C1212 = parameters[  6 ];
            C1213 = parameters[  7 ];
            C1222 = parameters[  8 ];
            C1223 = parameters[  9 ];
            C1233 = parameters[ 10 ];
            C1313 = parameters[ 11 ];
            C1322 = parameters[ 12 ];
            C1323 = parameters[ 13 ];
            C1333 = parameters[ 14 ];
            C2222 = parameters[ 15 ];
            C2223 = parameters[ 16 ];
            C2233 = parameters[ 17 ];
            C2323 = parameters[ 18 ];
            C2333 = parameters[ 19 ];
            C3333 = parameters[ 20 ];

        }
        else if ( parameters.size( ) == 9 ){

            C1111 = parameters[  0 ];
            C1122 = parameters[  1 ];
            C1133 = parameters[  2 ];
            C1212 = parameters[  3 ];
            C1313 = parameters[  4 ];
            C2222 = parameters[  5 ];
            C2233 = parameters[  6 ];
            C2323 = parameters[  7 ];
            C3333 = parameters[  8 ];

        }
        else if ( parameters.size( ) == 5 ){

            C1111 = parameters[  0 ];
            C1122 = parameters[  1 ];
            C1133 = parameters[  2 ];
            C1212 = 0.5 * ( C1111 - C1122 );
            C1313 = parameters[  3 ];
            C2222 = C1111;
            C2233 = C1133;
            C2323 = C1313;
            C3333 = parameters[  4 ];

        }
        else if ( parameters.size( ) == 3 ){

            C1111 = parameters[  0 ];
            C1122 = parameters[  1 ];
            C1133 = C1122;
            C1212 = parameters[  2 ];
            C1313 = C1212;
            C2222 = C1111;
            C2233 = C1122;
            C2323 = C1212;
            C3333 = C1111;

        }
        else if ( parameters.size( ) == 2 ){

            floatType lambda = parameters[ 0 ];
            floatType mu     = parameters[ 1 ];

            C1111 = lambda + 2 * mu;
            C1122 = lambda;
            C1133 = lambda;
            C1212 = mu;
            C1313 = C1212;
            C2222 = C1111;
            C2233 = lambda;
            C2323 = C1212;
            C3333 = C1111;

        }
        else{

            return new errorNode( __func__, "Requires 21 or 3 parameters. Parameters only defines " + std::to_string( parameters.size( ) ) );

        }

        stiffnessTensor = {
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

        return NULL;

    }

    errorOut rotateStiffnessTensor( const floatMatrix &directionCosines, const floatMatrix &stiffnessTensor,
                                    floatMatrix &rotatedStiffnessTensor ){
        /*!
         * Rotate the full 81 component stiffness tensor as
         *
         * \f$ C'_{ijkl} = R_{im} R_{jn} R_{ko} R_{lp} C_{mnop}
         *
         * where \f$C'_{ijkl}\f$ is the rotated stiffness tensor, \f$R_{ij}\f$ is the rotation matrix, and
         * \f$C_{mnop}\f$ is the original stiffness tensor.
         *
         * \param &directionCosines: The rotation matrix, \f$R_{ij}\f$
         * \param &stiffnessTensor: The stiffness tensor to rotate, \f$C_{mnop}\f$
         * \param &rotatedStiffnessTensor: The rotated stiffness tensor, \f$C'_{ijkl}\f$
         */

        rotatedStiffnessTensor = floatMatrix( spatialDimensions * spatialDimensions, floatVector( spatialDimensions * spatialDimensions, 0 ) );

        for ( unsigned int i=0; i<spatialDimensions; i++ ){
            for ( unsigned int j=0; j<spatialDimensions; j++ ){
                for ( unsigned int k=0; k<spatialDimensions; k++ ){
                    for ( unsigned int l=0; l<spatialDimensions; l++ ){
                        for ( unsigned int m=0; m<spatialDimensions; m++ ){
                            for ( unsigned int n=0; n<spatialDimensions; n++ ){
                                for ( unsigned int o=0; o<spatialDimensions; o++ ){
                                     for ( unsigned int p=0; p<spatialDimensions; p++ ){

            rotatedStiffnessTensor[ ( spatialDimensions * i ) + j ][ ( spatialDimensions * k ) + l ] +=
                directionCosines[ i ][ m ] * directionCosines[ j ][ n ] * directionCosines[ k ][ o ] * directionCosines[ l ][ p ]
                * stiffnessTensor[ ( spatialDimensions * m ) + n ][ ( spatialDimensions * o ) + p ];

                                     }
                                }
                            }
                        }
                    }
                }
            }
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
         * \param &parameters: The tensor components of the 9x9 stiffness tensor. Vector length determines the symmetry.
         *
         * - 81: A row-major vector representing the full 9x9 stiffness tensor directly.
         * - 21: fully anistropic \f$C_{1111}\f$, \f$C_{1112}\f$, \f$C_{1113}\f$, \f$C_{1122}\f$, \f$C_{1123}\f$,
         *   \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1213}\f$, \f$C_{1222}\f$, \f$C_{1223}\f$, \f$C_{1233}\f$,
         *   \f$C_{1313}\f$, \f$C_{1322}\f$, \f$C_{1323}\f$, \f$C_{1333}\f$, \f$C_{2222}\f$, \f$C_{2223}\f$,
         *   \f$C_{2233}\f$, \f$C_{2323}\f$, \f$C_{2333}\f$, \f$C_{3333}\f$
         * - 9: orthotropic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1313}\f$,
         *   \f$C_{2222}\f$, \f$C_{2323}\f$, \f$C_{3333}\f$
         * - 5: transversly isotropic or hexagonal \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1313}\f$, \f$C_{3333}\f$
         * - 3: cubic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1212}\f$
         * - 2: isotropic: lambda (\f$C_{1122}\f$), mu (\f$C_{1212}\f$)
         *
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
         * \param &parameters: The tensor components of the 9x9 stiffness tensor. Vector length determines the symmetry.
         *
         * - 81: A row-major vector representing the full 9x9 stiffness tensor directly.
         * - 21: fully anistropic \f$C_{1111}\f$, \f$C_{1112}\f$, \f$C_{1113}\f$, \f$C_{1122}\f$, \f$C_{1123}\f$,
         *   \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1213}\f$, \f$C_{1222}\f$, \f$C_{1223}\f$, \f$C_{1233}\f$,
         *   \f$C_{1313}\f$, \f$C_{1322}\f$, \f$C_{1323}\f$, \f$C_{1333}\f$, \f$C_{2222}\f$, \f$C_{2223}\f$,
         *   \f$C_{2233}\f$, \f$C_{2323}\f$, \f$C_{2333}\f$, \f$C_{3333}\f$
         * - 9: orthotropic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1313}\f$,
         *   \f$C_{2222}\f$, \f$C_{2323}\f$, \f$C_{3333}\f$
         * - 5: transversly isotropic or hexagonal \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1313}\f$, \f$C_{3333}\f$
         * - 3: cubic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1212}\f$
         * - 2: isotropic: lambda (\f$C_{1122}\f$), mu (\f$C_{1212}\f$)
         *
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
         * \param &parameters: The tensor components of the 9x9 stiffness tensor. Vector length determines the symmetry.
         *
         * - 81: A row-major vector representing the full 9x9 stiffness tensor directly.
         * - 21: fully anistropic \f$C_{1111}\f$, \f$C_{1112}\f$, \f$C_{1113}\f$, \f$C_{1122}\f$, \f$C_{1123}\f$,
         *   \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1213}\f$, \f$C_{1222}\f$, \f$C_{1223}\f$, \f$C_{1233}\f$,
         *   \f$C_{1313}\f$, \f$C_{1322}\f$, \f$C_{1323}\f$, \f$C_{1333}\f$, \f$C_{2222}\f$, \f$C_{2223}\f$,
         *   \f$C_{2233}\f$, \f$C_{2323}\f$, \f$C_{2333}\f$, \f$C_{3333}\f$
         * - 9: orthotropic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1313}\f$,
         *   \f$C_{2222}\f$, \f$C_{2323}\f$, \f$C_{3333}\f$
         * - 5: transversly isotropic or hexagonal \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1313}\f$, \f$C_{3333}\f$
         * - 3: cubic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1212}\f$
         * - 2: isotropic: lambda (\f$C_{1122}\f$), mu (\f$C_{1212}\f$)
         *
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
         * \param &parameters: The tensor components of the 9x9 stiffness tensor. Vector length determines the symmetry.
         *
         * - 81: A row-major vector representing the full 9x9 stiffness tensor directly.
         * - 21: fully anistropic \f$C_{1111}\f$, \f$C_{1112}\f$, \f$C_{1113}\f$, \f$C_{1122}\f$, \f$C_{1123}\f$,
         *   \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1213}\f$, \f$C_{1222}\f$, \f$C_{1223}\f$, \f$C_{1233}\f$,
         *   \f$C_{1313}\f$, \f$C_{1322}\f$, \f$C_{1323}\f$, \f$C_{1333}\f$, \f$C_{2222}\f$, \f$C_{2223}\f$,
         *   \f$C_{2233}\f$, \f$C_{2323}\f$, \f$C_{2333}\f$, \f$C_{3333}\f$
         * - 9: orthotropic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1212}\f$, \f$C_{1313}\f$,
         *   \f$C_{2222}\f$, \f$C_{2323}\f$, \f$C_{3333}\f$
         * - 5: transversly isotropic or hexagonal \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1133}\f$, \f$C_{1313}\f$, \f$C_{3333}\f$
         * - 3: cubic \f$C_{1111}\f$, \f$C_{1122}\f$, \f$C_{1212}\f$
         * - 2: isotropic: lambda (\f$C_{1122}\f$), mu (\f$C_{1212}\f$)
         *
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

}  // linearElasticity
}  // stressTools
