/**
  ******************************************************************************
  * \file tardigrade_mass_change_deformation.h
  ******************************************************************************
  * The header file for the deformation that results from mass change.
  ******************************************************************************
  */

#ifndef TARDIGRADE_MASS_CHANGE_DEFORMATION_H
#define TARDIGRADE_MASS_CHANGE_DEFORMATION_H

#define USE_EIGEN

#include<vector>
#include<array>

namespace tardigradeStressTools{

    namespace massChangeDeformation{

        constexpr int spatial_dimension = 3;

        typedef double floatType;
        typedef std::array< floatType, spatial_dimension > vector3d;
        typedef std::array< floatType, spatial_dimension > secondOrderTensor;

        //! Base class for the calculation of the deformation associated with a change in mass
        template< std::size_t num_params = 2 >
        class massChangeDeformationBase{

            public:

                massChangeDeformationBase( const floatType &dt,     const secondOrderTensor &At, const floatType &ct,     const floatType &ctp1,
                                           const floatType &rhot,   const floatType &rhotp1,     const floatType &gammat,
                                           const std::array< floatType, num_params > &parameters,
                                           const floatType &alpha=0.5, const floatType &tolr=1e-9, const floatType &tola=1e-9, const unsigned int &maxiter=20 );

            protected:

                const floatType _dt; //!< The change in time

                const secondOrderTensor _At; //!< The previous mass-deformation tensor

                const floatType _ct; //!< The previous density change rate

                const floatType _ctp1; //!< The current density change rate

                const floatType _rhot; //!< The previous density

                const floatType _rhotp1; //!< The current density

                const floatType _gammat; //!< The previous increment's rate multiplier

                const std::array< floatType, num_params > _parameters; //!< The model parameters

                const floatType &_alpha; //!< The integration parameter (0 for explicit, 1 for implicit)

                const floatType &_tolr; //!< The relative tolerance

                const floatType &_tola; //!< The absolute tolerance

                const unsigned int &_maxiter; //!< The maximuim number of iterations

                virtual void solveMassJacobian( );

                virtual void computeFlowDirection( );

                virtual void computeGammaRHS( );

                virtual void computeGammaLHS( );

                virtual void solveForGamma( );

            private:

                floatType _JAt;

                floatType _JAtp1;

                floatType _dJAtp1dC;

                floatType _dJAtp1dRho;

                secondOrderTensor _nt;

                secondOrderTensor _ntp1;

                floatType _gammatp1;

                floatType _gammaRHS;

                floatType _gammaLHS;

                secondOrderTensor _dGammaRHSdNtp1;

                floatType _dGammaRHSdJAtp1;

                secondOrderTensor _dGammadJAtp1;

                secondOrderTensor _dGammadNtp1;

        };

    }

}

#include "tardigrade_mass_change_deformation.cpp"

#endif
