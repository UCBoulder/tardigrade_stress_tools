/**
  ******************************************************************************
  * \file tardigrade_mass_change_deformation.cpp
  ******************************************************************************
  * The source file for the deformation that results from mass change.
  ******************************************************************************
  */

#include "tardigrade_mass_change_deformation.h"

#include "tardigrade_vector_tools.h"

namespace tardigradeStressTools{

    namespace massChangeDeformation{

        template< std::size_t num_params >
        massChangeDeformationBase<num_params>::massChangeDeformationBase( const floatType &dt,     const secondOrderTensor &At, const floatType &ct,     const floatType &ctp1,
                                                                          const floatType &rhot,   const floatType &rhotp1,
                                                                          const floatType &gammat, const std::array< floatType, num_params > &parameters,
                                                                          const floatType &alpha,  const floatType &tolr, const floatType &tola, const unsigned int &maxiter )
                                                                          : _dt( dt ), _At( At ), _ct( ct ), _ctp1( ctp1 ), _rhot( rhot ), _rhotp1( rhotp1 ),
                                                                            _gammat( gammat ),
                                                                            _parameters( parameters ), _alpha( alpha ), _tolr( tolr ), _tola( tola ), _maxiter( maxiter ){
            /*!
             * The base class for the mass change deformation family of material models where the 
             * evolution of the mass-change deformation is described by the form
             *
             * \f$ \ell_{ij}^{A} = \gamma n_{ij} \f$
             *
             * where $\bf{\ell}^A \f$ is the mass-change velocity gradient, \f$ \gamma \f$ is the 
             * rate multiplier, and \f$ n_{ij} \f$ is the direction tensor. This can be used to define
             * the ODE for the evolution of the mass-change deformation via
             *
             * \f$ \dot{A}_{iI} = \ell_{ij} A_{iI} \f$
             *
             * Expects to perform calculations in 3d
             *
             * \param &dt: The change in time (units: \f$\frac{m}{dv}\f$ )
             * \param &At: The previous value of the mass change deformation gradient (units: None, size=9)
             * \param &ct: The previous value of the mass density rate of change (units: \f$\frac{m}{dv t}\f$)
             * \param &ctp1: The current value of the mass density rate of change (units: \f$\frac{m}{dv t}\f$)
             * \param &rhot: The previous value of the mass density (units: \f$ \frac{m}{dv} \f$)
             * \param &rhotp1: The current value of the mass density (units: \f$ \frac{m}{dv} \f$)
             * \param &gammat: The previous rate multiplier for the evolution of the mass-change deformation (units: None)
             * \param &parameters: The parameter vector
             * \param &alpha: The integration parameter (0 for explicit 1 for implicit). Defaults to the second-order
             *    accurate value of 0.5
             */

            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > _At_map( _At.data( ) );

            _JAt = _At_map.determinant( );

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::solveMassJacobian( ){
            /*!
             * Solve for the updated mass Jacobian numerically using the midpoint rule
             */

            floatType rt = _ct / _rhot;

            floatType rtp1 = _ctp1 / _rhotp1;

            floatType num = 1. + _dt * ( 1. - _alpha ) * rt;

            floatType den = 1. - _dt * _alpha * rtp1;

            _JAtp1 = ( num / den ) * _JAt;

            _dJAtp1dC = num / ( den * den ) * _dt * _alpha / _rhotp1 * _JAt;

            _dJAtp1dRho = -num / ( den * den ) * _dt * _alpha * _ctp1 / ( _rhotp1 * _rhotp1 ) * _JAt;

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::computeFlowDirection( ){
            /*!
             * Compute the flow direction of the mass-change deformation
             */

            _nt = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

            _ntp1 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::computeGammaRHS( ){
            /*!
             * Compute the right-hand side for the iteration to solve for gammatp1
             */

            secondOrderTensor term1, term2;

            std::fill( std::begin( term1 ), std::end( term1 ), 0 );

            std::fill( std::begin( term2 ), std::end( term2 ), 0 );

            floatType factor = _JAt / _JAtp1;

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){ term1[ spatial_dimension * i + i ] = 1; term2[ spatial_dimension * i + i ] = factor; }

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){

                for ( unsigned int j = 0; j < spatial_dimension; j++ ){

                    term1[ spatial_dimension * i + j ] -= _alpha * _dt * _gammatp1 * _ntp1[ spatial_dimension * i + j ];
                    term2[ spatial_dimension * i + j ] += ( 1. - _alpha ) * _dt * _gammat * _nt[ spatial_dimension * i + j ] * factor;

                }

            }

            // Compute the determinant
            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > term1_map( term1.data( ) );

            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > term2_map( term2.data( ) );

            floatType _gammaRHS = term1_map.determinant( ) - term2_map.determinant( );

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::computeGammaLHS( ){
            /*!
             * Compute the left-hand side for the iteration to solve for gammatp1
             */

            secondOrderTensor term1, term2, invterm1;

            std::fill( std::begin( term1 ), std::end( term1 ), 0 );

            std::fill( std::begin( term2 ), std::end( term2 ), 0 );

            std::fill( std::begin( invterm1 ), std::end( invterm1 ), 0 );

            floatType factor = _JAt / _JAtp1;

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){ term1[ spatial_dimension * i + i ] = 1; term2[ spatial_dimension * i + i ] = factor; }

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){

                for ( unsigned int j = 0; j < spatial_dimension; j++ ){

                    term1[ spatial_dimension * i + j ] -= _alpha * _dt * _gammatp1 * _ntp1[ spatial_dimension * i + j ];
                    term2[ spatial_dimension * i + j ] += ( 1. - _alpha ) * _dt * _gammat * _nt[ spatial_dimension * i + j ] * factor;

                }

            }

            // Compute the inverse of term1
            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > term1_map( term1.data( ) );
            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > term2_map( term2.data( ) );
            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > invterm1_map( invterm1.data( ) );

            floatType detTerm1 = term1_map.determinant( );

            _dGammaRHSdJAtp1 = -term2_map.determinant( ) / _JAtp1;

            invterm1_map = term1_map.inverse( ).eval( );

            _gammaLHS = 0.;

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){
                for ( unsigned int j = 0; j < spatial_dimension; j++ ){
                    _gammaLHS += -detTerm1 * invterm1[ spatial_dimension * j + i ] * _alpha * _dt * _ntp1[ spatial_dimension * i + j ];
                    _dGammaRHSdNtp1[ spatial_dimension * i + j ] += -detTerm1 * invterm1[ spatial_dimension * j + i ] * _alpha * _dt * _gammatp1;
                }
            }

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::solveForGamma( ){
            /*!
             * Solve for the rate multiplier at the current timestep
             */

            // Set the value of the initial iterate
            _gammatp1 = 0;

            // Compute the initial value of the residual
            computeGammaRHS( );

            // Set the tolerance
            floatType tol = _tolr * std::fabs( _gammaRHS ) + _tola;

            unsigned int niter = 0;

            while ( ( niter < _maxiter ) && ( std::fabs( _gammaRHS ) > tol ) ){

                computeGammaLHS( );

                _gammatp1 += -_gammaRHS / _gammaLHS;

                computeGammaRHS( );

                niter++;

            }

            if ( std::fabs( _gammaRHS ) > tol ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::exception( "The solve for gamma did not converge" ) );

            }

            // Compute the derivatives
            computeGammaLHS( );

            _dGammadJAtp1 = -_dGammaRHSdJAtp1 / _gammaLHS;
            for ( auto v = std::begin( _dGammaRHSdNtp1 ); v != std::end( _dGammaRHSdNtp1 ); v++ ){
                _dGammadNtp1[ ( unsigned int )( v - std::begin( _dGammaRHSdNtp1 ) ) ] = -( *v ) / _gammaLHS;
            }

        }

        template< std::size_t num_params >
        void massChangeDeformationBase<num_params>::computeMassDeformation( ){
            /*!
             * Compute the deformation associated with the mass change
             */

            secondOrderTensor LHS, TERM1, RHS, invLHS;

            std::fill( std::begin( LHS ), std::end( LHS ), 0 );

            std::fill( std::begin( TERM1 ), std::end( TERM1 ), 0 );

            std::fill( std::begin( RHS ), std::end( RHS ), 0 );

            std::fill( std::begin( invLHS ), std::end( invLHS ), 0 );

            std::fill( std::begin( _Atp1 ), std::end( _Atp1 ), 0 );

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){ LHS[ spatial_dimension * i + i ] = 1.; TERM1[ spatial_dimension * i + i ] = 1.; }

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){

                for ( unsigned int j = 0; j < spatial_dimension; j++ ){

                    LHS[ spatial_dimension * i + j ] -= _alpha * _dt * _gammatp1 * _ntp1[ spatial_dimension * i + j ];

                    TERM1[ spatial_dimension * i + j ] += ( 1. - _alpha ) * _dt * _gammat * _nt[ spatial_dimension * i + j ];

                }

            }

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){

                for ( unsigned int j = 0; j < spatial_dimension; j++ ){

                    for ( unsigned int k = 0; k < spatial_dimension; k++ ){

                        RHS[ spatial_dimension * i + k ] += TERM1[ spatial_dimension * i + j ] * _At[ spatial_dimension * j + k ];

                    }

                }

            }

            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > LHS_map( LHS.data( ) );
            Eigen::Map< const Eigen::Matrix< floatType, spatial_dimension, spatial_dimension, Eigen::RowMajor > > invLHS_map( invLHS.data( ) );
            invLHS_map = LHS_map.inverse( ).eval( );

            for ( unsigned int i = 0; i < spatial_dimension; i++ ){

                for ( unsigned int j = 0; j < spatial_dimension; j++ ){

                    for ( unsigned int k = 0; k < spatial_dimension; k++ ){

                        _Atp1[ spatial_dimension * i + k ] += invLHS[ spatial_dimension * i + j ] * RHS[ spatial_dimension * j + k ];

                    }

                }

            }

        }

    }

}
