/**
  ******************************************************************************
  * \file tardigrade_stress_tools.cpp
  ******************************************************************************
  *  A collection of tools which implement and solve stress-strain relationships
  *  in such a way to enable more rapid development of constitutive models which
  *  have capabilities which may not be contained within a collection of
  *  of constitutive models.
  ******************************************************************************
  */

#include<tardigrade_stress_tools.h>

namespace tardigradeStressTools{

    template<
        unsigned int dim,
        class stress_iterator, typename stress_type
    >
    void calculateMeanStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end,
        stress_type &meanStress
    ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor
         * \param &stress_end: The stopping iterator of the stress tensor
         * \param &meanStress: The mean stress scalar
         */

        tardigradeVectorTools::rowMajorTrace<dim,dim>(
            stress_begin, stress_end, meanStress
        );

        meanStress /= 3.;

    }

    template<
        unsigned int dim,
        typename stress_type, class stress_iterator
    >
    stress_type calculateMeanStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end
    ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor
         * \param &stress_end: The stopping iterator of the stress tensor
         * \return meanStress: The mean stress scalar
         */

        stress_type meanStress;
        calculateMeanStress<dim>(
            stress_begin, stress_end, meanStress
        );

        return meanStress;

    }

    template<
        unsigned int dim,
        class stress_iterator, typename stress_type,
        class jacobian_iterator
    >
    void calculateMeanStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end,
        stress_type &meanStress,
        jacobian_iterator jacobian_begin, jacobian_iterator jacobian_end
    ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor
         * \param &stress_end: The stopping iterator of the stress tensor
         * \param &meanStress: The mean stress scalar
         * \param &jacobian_begin: The starting iterator of the row major jacobian w.r.t. the stress tensor
         * \param &jacobian_end: The stopping iterator of the row major jacobian w.r.t. the stress tensor
         */

        using jacobian_type = typename std::iterator_traits<jacobian_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( jacobian_end - jacobian_begin ) == dim * dim,
            "jacobian has a size of " + std::to_string( jacobian_end - jacobian_begin ) + " but should have a size of " + std::to_string( dim * dim )
        );

        calculateMeanStress<dim>(
            stress_begin, stress_end, meanStress
        );

        std::fill(
            jacobian_begin, jacobian_end, jacobian_type( )
        );

        for ( unsigned int i = 0; i < dim; ++i ){ *( jacobian_begin + dim * i + i ) = 1./3.; }

    }

    void calculateMeanStress( const floatVector &stress, floatType &meanStress ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateMeanStress< 3 >(
                std::begin( stress ), std::end( stress ), meanStress
            );
        }
        else if ( dim == 2 ){
            calculateMeanStress< 2 >(
                std::begin( stress ), std::end( stress ), meanStress
            );
        }
        else if ( dim == 1 ){
            calculateMeanStress< 1 >(
                std::begin( stress ), std::end( stress ), meanStress
            );
        }

        return;
    }

    floatType calculateMeanStress( const floatVector &stress ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress: The stress tensor
         * \return meanStress: The mean stress scalar
         */
        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            return calculateMeanStress< 3, floatType >(
                std::begin( stress ), std::end( stress )
            );
        }
        else if ( dim == 2 ){
            return calculateMeanStress< 2, floatType >(
                std::begin( stress ), std::end( stress )
            );
        }
        else if ( dim == 1 ){
            return calculateMeanStress< 1, floatType >(
                std::begin( stress ), std::end( stress )
            );
        }

        return 0; //Default value. You never want to get here!

    }

    void calculateMeanStress( const floatMatrix &stress, floatType &meanStress ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress: The stress tensor
         * \param &meanStress: The mean stress scalar
         */

        const unsigned int dim = stress.size( );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        floatVector _stress = tardigradeVectorTools::appendVectors( stress );

        if ( dim == 3 ){
            return calculateMeanStress< 3 >(
                std::begin( _stress ), std::end( _stress ), meanStress
            );
        }
        else if ( dim == 2 ){
            return calculateMeanStress< 2 >(
                std::begin( _stress ), std::end( _stress ), meanStress
            );
        }
        else if ( dim == 1 ){
            return calculateMeanStress< 1 >(
                std::begin( _stress ), std::end( _stress ), meanStress
            );
        }
    }

    floatType calculateMeanStress( const floatMatrix &stress ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in matrix format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress: The stress tensor
         * \return meanStress: The mean stress scalar
         */

        floatType meanStress;
        TARDIGRADE_ERROR_TOOLS_CATCH( calculateMeanStress( stress, meanStress ) )

        return meanStress;
    }

    void calculateMeanStress( const floatVector &stress, floatType &meanStress, floatVector &jacobian ){
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress: The row major stress tensor
         * \param &meanStress: The scalar mean stress
         * \param &jacobian: The row major jacobian w.r.t. the stress tensor
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        jacobian = floatVector( dim * dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateMeanStress< 3 >(
                std::begin( stress ),   std::end( stress ), meanStress,
                std::begin( jacobian ), std::end( jacobian )
            );
        }
        else if ( dim == 2 ){
            calculateMeanStress< 2 >(
                std::begin( stress ),   std::end( stress ), meanStress,
                std::begin( jacobian ), std::end( jacobian )
            );
        }
        else if ( dim == 1 ){
            calculateMeanStress< 1 >(
                std::begin( stress ),   std::end( stress ), meanStress,
                std::begin( jacobian ), std::end( jacobian )
            );
        }

        return;

    }

    template<
        unsigned int dim,
        class stress_iterator, class deviatoric_iterator
    >
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(
        const stress_iterator &stress_begin,  const stress_iterator &stress_end,
        deviatoric_iterator deviatoric_begin, deviatoric_iterator deviatoric_end
    ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor in row major format
         * \param &stress_end: The stopping iterator of the stress tensor in row major format
         * \param deviatoric_begin: The starting iterator of the deviatoric stress tensor in row major format
         * \param deviatoric_end: The stopping iterator of the deviatoric stress tensor in row major format
         */

        using stress_type = typename std::iterator_traits<stress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( stress_end - stress_begin ) == dim * dim,
            "The stress has a size of " + std::to_string( ( unsigned int )( stress_end - stress_begin ) ) + " but must have a size of " + std::to_string( dim * dim )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( deviatoric_end - deviatoric_begin ) == dim * dim,
            "The deviatoric stress has a size of " + std::to_string( ( unsigned int )( deviatoric_end - deviatoric_begin ) ) + " but must have a size of " + std::to_string( dim * dim )
        );

        stress_type meanStress;

        calculateMeanStress< dim >(
            stress_begin, stress_end, meanStress
        );

        std::copy(
            stress_begin, stress_end, deviatoric_begin
        );

        for ( unsigned int i = 0; i < dim; ++i ){ *( deviatoric_begin + dim * i + i ) -= meanStress; }

        return;

    }

    template<
        unsigned int dim,
        class stress_iterator, class deviatoric_iterator,
        class jacobian_iterator
    >
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(
        const stress_iterator &stress_begin,  const stress_iterator &stress_end,
        deviatoric_iterator deviatoric_begin, deviatoric_iterator deviatoric_end,
        jacobian_iterator jacobian_begin,     jacobian_iterator jacobian_end
    ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor in row major format
         * \param &stress_end: The stopping iterator of the stress tensor in row major format
         * \param deviatoric_begin: The starting iterator of the deviatoric stress tensor in row major format
         * \param deviatoric_end: The stopping iterator of the deviatoric stress tensor in row major format
         * \param jacobian_begin: The starting iterator of the row-major Jacobian
         * \param jacobian_end: The stopping iterator of the row-major Jacobian
         */

        using stress_type = typename std::iterator_traits<stress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( jacobian_end - jacobian_begin ) == dim * dim * dim * dim,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( jacobian_end - jacobian_begin ) ) + " but must have a size of " + std::to_string( dim * dim * dim * dim )
        );

        calculateDeviatoricStress< dim >(
            stress_begin, stress_end, deviatoric_begin, deviatoric_end
        );

        std::fill(
            jacobian_begin, jacobian_end, stress_type( )
        );

        for ( unsigned int i = 0; i < dim * dim; ++i ){

            *( jacobian_begin + dim * dim * i + i ) += 1;


        }

        for ( unsigned int i = 0; i < dim; ++i ){

            for ( unsigned int j = 0; j < dim; ++j ){

                *( jacobian_begin + dim * dim * dim * i + dim * dim * i + dim * j + j ) -= 1./3;

            }

        }

        return;

    }

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &deviatoric: The deviatoric stress tensor in row major format
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        deviatoric = floatVector( dim * dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateDeviatoricStress< 3 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric )
            );
        }
        else if ( dim == 2 ){
            calculateDeviatoricStress< 2 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric )
            );
        }
        else if ( dim == 1 ){
            calculateDeviatoricStress< 1 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric )
            );
        }

        return;

    }

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric, floatVector &jacobian ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * Also return the jacobian
         *
         * \f$\frac{ \partial \sigma_{ ij }^{ deviatoric } }{ \partial \sigma_{ kl } } = \delta_{ ik }\delta{ jl } - \frac{ \partial
         * \bar{ \sigma } }{ \partial \sigma_{ kl } } \delta_{ ij }\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &deviatoric: The deviatoric stress tensor in row major format
         * \param &jacobian: The jacobian of the deviatoric stress tensor w.r.t. the stress in row-major format
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        deviatoric = floatVector( dim * dim, 0 );
        jacobian   = floatVector( dim * dim * dim * dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateDeviatoricStress< 3 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric ),
                std::begin( jacobian ),   std::end( jacobian )
            );
        }
        else if ( dim == 2 ){
            calculateDeviatoricStress< 2 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric ),
                std::begin( jacobian ),   std::end( jacobian )
            );
        }
        else if ( dim == 1 ){
            calculateDeviatoricStress< 1 >(
                std::begin( stress ),     std::end( stress ),
                std::begin( deviatoric ), std::end( deviatoric ),
                std::begin( jacobian ),   std::end( jacobian )
            );
        }

        return;

    }

    void calculateDeviatoricStress( const floatVector &stress, floatVector &deviatoric, floatMatrix &jacobian ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * Also return the jacobian
         *
         * \f$\frac{ \partial \sigma_{ ij }^{ deviatoric } }{ \partial \sigma_{ kl } } = \delta_{ ik }\delta{ jl } - \frac{ \partial
         * \bar{ \sigma } }{ \partial \sigma_{ kl } } \delta_{ ij }\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &deviatoric: The deviatoric stress tensor in row major format
         * \param &jacobian: The jacobian of the deviatoric stress tensor w.r.t. the stress.
         */

        floatVector _jacobian;

        TARDIGRADE_ERROR_TOOLS_CATCH( calculateDeviatoricStress( stress, deviatoric, _jacobian ) );

        jacobian = tardigradeVectorTools::inflate( _jacobian, stress.size( ), stress.size( ) ); 

        return;

    }

    floatVector calculateDeviatoricStress( const floatVector &stress ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress: The stress tensor in row major format
         * \return deviatoric: The deviatoric stress tensor in row major format
         */

        floatVector deviatoric( stress.size( ) );
        TARDIGRADE_ERROR_TOOLS_CATCH( calculateDeviatoricStress( stress, deviatoric ) )

        return deviatoric;
    }

    floatVector calculateDeviatoricStress( const floatVector &stress, floatMatrix &jacobian ){
        /*!
         * Compute the deviatoric stress tensor from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * Also return the jacobian
         *
         * \f$\frac{ \partial \sigma_{ ij }^{ deviatoric } }{ \partial \sigma_{ kl } } = \delta_{ ik }\delta{ jl } - \frac{ \partial
         * \bar{ \sigma } }{ \partial \sigma_{ kl } } \delta_{ ij }\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &jacobian: The jacobian of the deviatoric stress tensor w.r.t. the stress.
         * \return deviatoric: The deviatoric part of the stress tensor in row major format
         */

        floatVector deviatoric, _jacobian;

        TARDIGRADE_ERROR_TOOLS_CATCH( calculateDeviatoricStress( stress, deviatoric, _jacobian ) );

        jacobian = tardigradeVectorTools::inflate( _jacobian, stress.size( ), stress.size( ) ); 

        return deviatoric;

    }

    template<
        unsigned int dim,
        class stress_iterator, typename vonMises_type
    >
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end,
        vonMises_type &vonMises
    ){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ vonMises } = \sqrt{ \frac{ 3 }{ 2 }*\sigma^{ deviatoric }\sigma^{ deviatoric } }\f$
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress_begin: The starting iterator of the row major stress tensor
         * \param &stress_end: The stopping iterator of the row major stress tensor
         * \param &vonMises: The von-Mises stress
         */

        using stress_type = typename std::iterator_traits<stress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( stress_end - stress_begin ) == dim * dim,
            "The stress tensor has a size of " + std::to_string( ( unsigned int )( stress_end - stress_begin ) ) + " but should have a size of " + std::to_string( dim * dim )
        );

        std::array< stress_type, dim * dim > deviatoric;
        calculateDeviatoricStress<dim>(
            stress_begin, stress_end,
            std::begin( deviatoric ), std::end( deviatoric )
        );

        vonMises = std::sqrt(
            1.5 * std::inner_product( std::begin( deviatoric ), std::end( deviatoric ), std::begin( deviatoric ), stress_type( ) )
        );

    }

    template<
        unsigned int dim,
        class stress_iterator, typename vonMises_type,
        class jacobian_iterator
    >
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end,
        vonMises_type &vonMises,
        jacobian_iterator jacobian_begin, jacobian_iterator jacobian_end
    ){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ vonMises } = \sqrt{ \frac{ 3 }{ 2 }*\sigma^{ deviatoric }\sigma^{ deviatoric } }\f$
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress_begin: The starting iterator of the row major stress tensor
         * \param &stress_end: The stopping iterator of the row major stress tensor
         * \param &vonMises: The von-Mises stress
         * \param jacobian_begin: The starting iterator of the Jacobian
         * \param jacobian_end: The stopping iterator of the Jacobian
         */

        using stress_type = typename std::iterator_traits<stress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( stress_end - stress_begin ) == dim * dim,
            "The stress tensor has a size of " + std::to_string( ( unsigned int )( stress_end - stress_begin ) ) + " but should have a size of " + std::to_string( dim * dim )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( stress_end - stress_begin ) == dim * dim,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( jacobian_end - jacobian_begin ) ) + " but should have a size of " + std::to_string( dim * dim )
        );

        std::array< stress_type, dim * dim > deviatoric;
        calculateDeviatoricStress<dim>(
            stress_begin, stress_end,
            std::begin( deviatoric ), std::end( deviatoric )
        );

        vonMises = std::sqrt(
            1.5 * std::inner_product( std::begin( deviatoric ), std::end( deviatoric ), std::begin( deviatoric ), stress_type( ) )
        );

        tardigradeConstitutiveTools::computeUnitNormal(
            std::begin( deviatoric ), std::end( deviatoric ),
            jacobian_begin,           jacobian_end
        );

        std::transform(
            jacobian_begin, jacobian_end, jacobian_begin,
             std::bind( std::multiplies< >( ), std::placeholders::_1, std::sqrt( 1.5 ) ) 
        );

    }

    void calculateVonMisesStress( const floatVector &stress, floatType &vonMises ){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ vonMises } = \sqrt{ \frac{ 3 }{ 2 }*\sigma^{ deviatoric }\sigma^{ deviatoric } }\f$
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress: The row major stress tensor
         * \param &vonMises: The von-Mises stress
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateVonMisesStress< 3 >(
                std::begin( stress ),     std::end( stress ),
                vonMises
            );
        }
        else if ( dim == 2 ){
            calculateVonMisesStress< 2 >(
                std::begin( stress ),     std::end( stress ),
                vonMises
            );
        }
        else if ( dim == 1 ){
            calculateVonMisesStress< 1 >(
                std::begin( stress ),     std::end( stress ),
                vonMises
            );
        }

        return;

    }

    floatType calculateVonMisesStress( const floatVector &stress ){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ vonMises } = \sqrt{ \frac{ 3 }{ 2 }*\sigma^{ deviatoric }\sigma^{ deviatoric } }\f$
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress: The row major stress tensor
         * \return vonMises: The von-Mises stress
         */

        floatType vonMises = 0.;
        TARDIGRADE_ERROR_TOOLS_CATCH( calculateVonMisesStress( stress, vonMises ) )

        return vonMises;
    }

    void calculateVonMisesStress( const floatVector &stress, floatType &vonMises, floatVector &jacobian ){
        /*!
         * Compute the von Mises stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$\sigma^{ vonMises } = \sqrt{ \frac{ 3 }{ 2 }*\sigma^{ deviatoric }\sigma^{ deviatoric } }\f$
         *
         * \f$\sigma^{ deviatoric } = \sigma - \sigma^{ mean }I\f$
         *
         * \param &stress: The row major stress tensor
         * \param &vonMises: The von mises stress
         * \param &jacobian: The row major mean stress jacobian tensor w.r.t. the stress tensor
         */

        const unsigned int dim = std::pow( stress.size( ), 0.5 );

        jacobian = floatVector( dim * dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( dim == 3 ) || ( dim == 2 ) || ( dim == 1 ),
            "The dimension of the stress is " + std::to_string( dim ) + " but must be 1, 2, or 3"
        );

        if ( dim == 3 ){
            calculateVonMisesStress< 3 >(
                std::begin( stress ),   std::end( stress ),
                vonMises,
                std::begin( jacobian ), std::end( jacobian )
            );
        }
        else if ( dim == 2 ){
            calculateVonMisesStress< 2 >(
                std::begin( stress ),     std::end( stress ),
                vonMises,
                std::begin( jacobian ), std::end( jacobian )
            );
        }
        else if ( dim == 1 ){
            calculateVonMisesStress< 1 >(
                std::begin( stress ),     std::end( stress ),
                vonMises,
                std::begin( jacobian ), std::end( jacobian )
            );
        }

        return;

    }

//    template<
//        typename vonMises_type, typename meanStress_type,
//        class dpParam_iterator, typename dpYield_type
//    >
//    void druckerPragerSurface(
//        const vonMises_type &vonMises, const meanStress_type &meanStress,
//        const dpParam_iterator &dpParam_begin, const dpParam_iterator &dpParam_end,
//        dpYield_type &dpYield
//    ){
//        /*!
//         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
//         *
//         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
//         *
//         * TODO: find the common name for which material parameter, if a common
//         * name exists to distinguish between the two DP parameters.
//         *
//         * \param &vonMises: The von Mises stress
//         * \param &meanStress: The mean Stress
//         * \param &dpParam_begin: The starting iterator of the Drucker-Prager material parameters
//         * \param &dpParam_end: The stopping iterator of the Drucker-Prager material parameters
//         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
//         */
//
//        TARDIGRADE_ERROR_TOOLS_CHECK(
//            ( unsigned int )( dpParam_end - dpParam_begin ) == 2,
//            "The Drucker-Prager parameter vector has a size of " + std::to_string( ( unsigned int )( dpParam_end - dpParam_begin ) ) + " but must have a size of 2"
//        );
//
//        druckerPragerSurface(
//            vonMises, meanStress,
//            *( dpParam_begin + 0 ), *( dpParam_begin + 1 ),
//            dpYield
//        );
//
//    }
//
//    template<
//        typename vonMises_type, typename meanStress_type,
//        typename A_type, typename B_type,
//        typename dpYield_type
//    >
//    dpYield_type druckerPragerSurface(
//        const vonMises_type &vonMises, const meanStress_type &meanStress,
//        const A_type &A, const B_type &B, dpYield_type &dpYield
//    ){
//
//    }
//
//    template<
//        typename vonMises_type, typename meanStress_type,
//        class dpParam_iterator, typename dpYield_type
//    >
//    dpYield_type druckerPragerSurface(
//        const vonMises_type &vonMises, const meanStress_type &meanStress,
//        const dpParam_iterator &dpParam_begin, const dpParam_iterator &dpParam_end,
//        dpYield_type &dpYield
//    ){
//
//    }

//    void druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B, floatType &dpYield ){
//        /*!
//         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
//         *
//         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
//         *
//         * TODO: find the common name for which material parameter, if a common
//         * name exists to distinguish between the two DP parameters.
//         *
//         * \param &vonMises: The von Mises stress
//         * \param &meanStress: The mean Stress
//         * \param &A: The first Drucker-Prager material parameter
//         * \param &B: The second Drucker-Prager material parameter
//         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
//         */
//
//        dpYield = vonMises + A*meanStress - B;
//
//        return;
//    }

    void druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam, floatType &dpYield ){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
         *
         * \f$f = \sigma^{ vonMises } + dpParam[ 0 ]*\sigma^{ mean } - dpParam[ 1 ]\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        //Check Drucker-Prager parameter vector length
        TARDIGRADE_ERROR_TOOLS_CHECK( dpParam.size( ) == 2, "Two parameters are required for the Drucker-Prager surface." );

        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( vonMises, meanStress, dpParam[ 0 ], dpParam[ 1 ], dpYield ) )

        return;
    }

    floatType druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatType &A, const floatType &B ){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \return dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0;

        //Calculate DP yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( vonMises, meanStress, A, B, dpYield ) );

        return dpYield;
    }

    floatType druckerPragerSurface( const floatType &vonMises, const floatType &meanStress, const floatVector &dpParam ){
        /*!
         * Compute the Drucker-Prager yield criterion from the von Mises and mean stress
         *
         * \f$f = \sigma^{ vonMises } + dpParam[ 0 ]*\sigma^{ mean } - dpParam[ 1 ]\f$
         *
         * \param &vonMises: The von Mises stress
         * \param &meanStress: The mean Stress
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \return dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0;
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( vonMises, meanStress, dpParam, dpYield ) );

        return dpYield;
    }

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
        calculateVonMisesStress( stress, vonMises );
        calculateMeanStress( stress, meanStress );

        //Calculate DP yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( vonMises, meanStress, A, B, dpYield ) )

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * \param &stress: The stress tensor in row major format
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        //Check Drucker-Prager parameter vector length
        TARDIGRADE_ERROR_TOOLS_CHECK( dpParam.size( ) == 2, "Two parameters are required for the Drucker-Prager surface." )

        //Calculate DP yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, dpParam[ 0 ], dpParam[ 1 ], dpYield ) )

        return;
    }

    floatType druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The row major stress tensor
         * \param &A: The first Drucker-Prager material parameter
         * \param &B: The second Drucker-Prager material parameter
         * \return dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, A, B, dpYield ) );

        return dpYield;
    }

    floatType druckerPragerSurface( const floatVector &stress, const floatVector &dpParam ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * \param &stress: The row major stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \return dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        floatType dpYield = 0.;

        //Calculate DP yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, dpParam, dpYield ) );

        return dpYield;
    }

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
        unsigned int length = stress.size( );
        TARDIGRADE_ERROR_TOOLS_CHECK( length == jacobian.size( ), "The stress tensor and jacobian tensor sizes must match." );

        //Calculate von Mises and mean stresses with jacobians
        floatType vonMises, meanStress;
        vonMises = meanStress = 0.;
        floatVector vonMisesJacobian( jacobian.size( ), 0. );
        floatVector meanStressJacobian( jacobian.size( ), 0. );
        calculateVonMisesStress( stress, vonMises, vonMisesJacobian );
        calculateMeanStress( stress, meanStress, meanStressJacobian );

        //Calculate the Drucker-Prager yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, A, B, dpYield ) )

        //Calculate the Drucker-Prager jacobian
        jacobian = vonMisesJacobian + A * meanStressJacobian;

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         */

        //Check Drucker-Prager parameter vector length
        TARDIGRADE_ERROR_TOOLS_CHECK( dpParam.size( ) == 2, "Two parameters are required for the Drucker-Prager surface." );

        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, dpParam[ 0 ], dpParam[ 1 ], dpYield, jacobian ) )

        return;
    }

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &djacobiandstress, const floatType tol ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
         * \param tol: The tolerance on the von-mises stress to prevent nan when the
         *     incoming stress is very small (defaults to 1e-9)
         */

        //Check vector lengths
        unsigned int length = stress.size( );
        TARDIGRADE_ERROR_TOOLS_CHECK( length == jacobian.size( ), "The stress tensor and jacobian tensor sizes must match." );

        //Calculate von Mises and mean stresses with jacobians
        floatType vonMises, meanStress;
        vonMises = meanStress = 0.;
        floatVector vonMisesJacobian( jacobian.size( ), 0. );
        floatVector meanStressJacobian( jacobian.size( ), 0. );
        calculateVonMisesStress( stress, vonMises, vonMisesJacobian );
        calculateMeanStress( stress, meanStress, meanStressJacobian );

        //Calculate the Drucker-Prager yield criterion
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, A, B, dpYield ) );

        //Calculate the Drucker-Prager jacobian
        jacobian = vonMisesJacobian + A * meanStressJacobian;

        //Compute the gradient of the jacobian w.r.t. the stress
        floatVector eye( stress.size( ), 0 );
        tardigradeVectorTools::eye<floatType>( eye );

        floatVector deviatoric = calculateDeviatoricStress( stress );

        const unsigned int rows = jacobian.size( );
        const unsigned int cols = stress.size( );
        djacobiandstress = floatVector( jacobian.size( ) * stress.size( ), 0 );

        if ( tardigradeVectorTools::l2norm<floatType>(std::begin( vonMisesJacobian ), std::end( vonMisesJacobian ) ) < tol ){
            return;
        }
        for ( unsigned int ij = 0; ij < rows; ++ij ){

            djacobiandstress[ cols * ij + ij ] += ( 3 / ( 2 * ( vonMises + tol ) ) );

            for ( unsigned int ab = 0; ab < cols; ++ab ){

                djacobiandstress[ cols * ij + ab ] -= ( 3 / ( 2 * ( vonMises + tol ) ) ) * ( eye[ ij ] * meanStressJacobian[ ab ] + deviatoric[ ij ] * vonMisesJacobian[ ab ] / ( vonMises + tol ) );

            }

        }

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress, const floatType tol ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
         * \param tol: The tolerance on the von-mises stress to prevent nan when the
         *     incoming stress is very small (defaults to 1e-9)
         */

        floatVector _djacobiandstress;
        TARDIGRADE_ERROR_TOOLS_CATCH(
            druckerPragerSurfaceFlatJ( stress, A, B, dpYield, jacobian, _djacobiandstress, tol );
        );

        djacobiandstress = tardigradeVectorTools::inflate(
            _djacobiandstress, jacobian.size( ), stress.size( )
        );

    }

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &djacobiandstress ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &djacobiandstress: The gradient of the jacobian w.r.t. the stress
         */

        //Check Drucker-Prager parameter vector length
        TARDIGRADE_ERROR_TOOLS_CHECK( dpParam.size( ) == 2, "Two parameters are required for the Drucker-Prager surface." )

        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, dpParam[ 0 ], dpParam[ 1 ], dpYield, jacobian, djacobiandstress ) )

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatMatrix &djacobiandstress ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * TODO: find the common name for which material parameter, if a common
         * name exists to distinguish between the two DP parameters.
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &djacobiandstress: The gradient of the jacobian w.r.t. the stress
         */

        //Check Drucker-Prager parameter vector length
        TARDIGRADE_ERROR_TOOLS_CHECK( dpParam.size( ) == 2, "Two parameters are required for the Drucker-Prager surface." )

        floatVector _djacobiandstress;
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, dpParam[ 0 ], dpParam[ 1 ], dpYield, jacobian, _djacobiandstress ) )

        djacobiandstress = tardigradeVectorTools::inflate( _djacobiandstress, jacobian.size( ), stress.size( ) );

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, const floatType tol ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
         * \param &tol: The tolerance of the unit direction norm to prevent nan
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, A, B, dpYield, jacobian ) )

        //Calculate the Drucker-Prager unit normal flow direction as the normalized jacobian
        unitDirection = floatVector( jacobian.size( ), 0 );
        auto norm = tardigradeVectorTools::l2norm<floatType>( std::begin( jacobian ), std::end( jacobian ) );
        std::transform(
            std::begin( jacobian ),
            std::end( jacobian ),
            std::begin( unitDirection ),
            std::bind(
                std::multiplies<floatType>(),
                std::placeholders::_1, 1./(norm+tol)
            )
        );

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, const floatType tol ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * \param &tol: The tolerance of the unit direction norm to prevent nan
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurface( stress, dpParam, dpYield, jacobian ) )

        //Calculate the Drucker-Prager unit normal flow direction as the normalized jacobian
        auto norm = tardigradeVectorTools::l2norm<floatType>( std::begin( jacobian ), std::end( jacobian ) );
        std::transform(
            std::begin( jacobian ),
            std::end( jacobian ),
            std::begin( unitDirection ),
            std::bind(
                std::multiplies<floatType>(),
                std::placeholders::_1, 1./(norm+tol)
            )
        );

        return;
    }

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatVector &unitDirectionJacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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
        floatVector djacobiandstress;
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, A, B, dpYield, jacobian, djacobiandstress ) )

        //Compute the unit normal flow direction and the jacobian of the unit normal flow direction
        //w.r.t. stress
        floatVector duDdjacobian;
        tardigradeConstitutiveTools::computeUnitNormal( jacobian, unitDirection, duDdjacobian );

        const unsigned int rows = unitDirection.size( );
        const unsigned int cols = stress.size( );
        unitDirectionJacobian = floatVector( rows * cols, 0 );

        if ( tardigradeVectorTools::l2norm<floatType>( std::begin( unitDirection ), std::end( unitDirection ) ) < 0.5 ){
            return;
        }

        const unsigned int jsize = jacobian.size( );
        for ( unsigned int i = 0; i < rows; ++i ){
            for ( unsigned int j = 0; j < jsize; ++j ){
                for ( unsigned int k = 0; k < cols; ++k ){
                    unitDirectionJacobian[ cols * i + k ] += duDdjacobian[ jsize * i + j ] * djacobiandstress[ cols * j + k ];
                }
            }
        }

        return;
    }


    void druckerPragerSurface( const floatVector &stress, const floatType &A, const floatType &B, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
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

        floatVector _unitDirectionJacobian;

        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, A, B, dpYield, jacobian, unitDirection, _unitDirectionJacobian ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionJacobian = tardigradeVectorTools::inflate( _unitDirectionJacobian, unitDirection.size( ), stress.size( ) ) );

        return;
    }

    void druckerPragerSurfaceFlatJ( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatVector &unitDirectionJacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * \param &unitDirectionJacobian: The jacobian of the unit direction w.r.t. the stress
         */

        //Calculate the Drucker-Prager yield criterion and jacobian
        floatVector djacobiandstress;
        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, dpParam, dpYield, jacobian, djacobiandstress ) )

        //Compute the unit normal flow direction and the jacobian of the unit normal flow direction
        //w.r.t. stress
        floatVector duDdjacobian;
        tardigradeConstitutiveTools::computeUnitNormal( jacobian, unitDirection, duDdjacobian );

        const unsigned int rows = unitDirection.size( );
        const unsigned int cols = stress.size( );
        unitDirectionJacobian = floatVector( rows * cols, 0 );
        if ( tardigradeVectorTools::l2norm<floatType>( std::begin( unitDirection ), std::end( unitDirection ) ) < 0.5 ){
            return;
        }

        const unsigned int jsize = jacobian.size( );
        for ( unsigned int i = 0; i < rows; ++i ){
            for ( unsigned int j = 0; j < jsize; ++j ){
                for ( unsigned int k = 0; k < cols; ++k ){
                    unitDirectionJacobian[ cols * i + k ] += duDdjacobian[ jsize * i + j ] * djacobiandstress[ cols * j + k ];
                }
            }
        }

        return;
    }

    void druckerPragerSurface( const floatVector &stress, const floatVector &dpParam, floatType &dpYield, floatVector &jacobian, floatVector &unitDirection, floatMatrix &unitDirectionJacobian ){
        /*!
         * Compute the Drucker-Prager yield criterion from a 2nd rank stress tensor stored in row major format
         *
         * \f$f = \sigma^{ vonMises } + A*\sigma^{ mean } - B\f$
         *
         * \param &stress: The stress tensor
         * \param &dpParam: The two Drucker-Prager material parameters in a vector { A, B }
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         * \param &jacobian: The row major jacobian tensor w.r.t. the stress
         * \param &unitDirection: The normalized row major jacobian tensor w.r.t. the stress
         * \param &unitDirectionJacobian: The jacobian of the unit direction w.r.t. the stress
         */


        floatVector _unitDirectionJacobian;

        TARDIGRADE_ERROR_TOOLS_CATCH( druckerPragerSurfaceFlatJ( stress, dpParam, dpYield, jacobian, unitDirection, _unitDirectionJacobian ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionJacobian = tardigradeVectorTools::inflate( _unitDirectionJacobian, unitDirection.size( ), stress.size( ) ) );

        return;
    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &dStress, floatVector &stress, floatVector &currentStateVariables ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStrain: The previous value of strain
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &dStress: The computed change in stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &stress: The computed stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         */

        floatType dt = currentTime - previousTime;

        //Check the material parameters
        TARDIGRADE_ERROR_TOOLS_CHECK( materialParameters.size( ) != 0, "At least one material parameter must be provided." );

        //Set the number of Prony-Series terms
        TARDIGRADE_ERROR_TOOLS_CHECK( ( materialParameters.size( ) - 1 ) % 2 == 0, "An equal number of taus and Gs must be provided." );

        unsigned int nTerms = ( materialParameters.size( ) - 1 )/2;

        //Set the dimension of the strain
        unsigned int dim = currentStrain.size( );
        //Check the state variables
        TARDIGRADE_ERROR_TOOLS_CHECK( previousStateVariables.size( ) == nTerms*dim, "The number of previous state variables is not consistent with the strain size." );

        //Compute the infinite stress
        floatVector dStrain = currentStrain - previousStrain;
        dStress = materialParameters[ 0 ]*dStrain;
        stress = materialParameters[ 0 ]*currentStrain;

        //Set the initial value of the factor
        floatType factor;
        floatType taui;
        floatType Gi;
        floatVector Xip( dim, 0 ), Xic( dim, 0 );
        currentStateVariables.resize( 0 );

        std::vector< unsigned int > indices( dim, 0 );

        floatVector dXi( dim, 0 );

        for ( unsigned int i=1; i<nTerms+1; i++ ){
            //Get the parameters
            taui = materialParameters[ i ];
            Gi = materialParameters[ i+nTerms ];

            //Get the factor
            factor = dt/( taui + dt*( 1 - alpha )*currentRateModifier );

            //Set the indices of the previous values of the state variables
            for ( unsigned int j=dim*( i-1 ), k=0; j<dim*i; j++, k++ ){
                indices[ k ] = j;
            }

            //Get the previous values of the state variables
            tardigradeVectorTools::getValuesByIndex( previousStateVariables, indices, Xip );

            //Compute the new state-variable values
            floatVector dxi = factor * ( ( 1 - alpha ) * ( currentStrain - Xip ) * currentRateModifier
                                           + alpha * ( previousStrain - Xip ) * previousRateModifier );

            Xic = Xip + dxi;

            //Add the contribution to the stress
            dStress += Gi*( dStrain - dxi );
            stress += Gi*( currentStrain - Xic );

            //Save the new value of the state variable
            currentStateVariables = tardigradeVectorTools::appendVectors( { currentStateVariables, Xic } );
        }

        return;
    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStrain: The previous value of strain
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &stress: The computed stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         */

        floatVector dStress;

        return linearViscoelasticity( currentTime, currentStrain, previousTime, previousStrain,
                                      currentRateModifier, previousRateModifier,
                                      previousStateVariables, materialParameters,
                                      alpha, dStress, stress, currentStateVariables );
    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &dStress, floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain,
                                   floatVector &dstressdrateModifier ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &previousStrain: The previous value of strain
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &dStress: The computed change in stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &stress: The computed stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         * \param &dstressdstrain: The derivative of the stress w.r.t. the strain
         * \param &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( linearViscoelasticity( currentTime, currentStrain, previousTime, previousStrain,
                                                             currentRateModifier, previousRateModifier, previousStateVariables,
                                                             materialParameters, alpha, dStress, stress, currentStateVariables ) )

        //Compute the change in time
        float dt = currentTime - previousTime;

        //Compute the derivative of the stress w.r.t. the strain
        //and the rate modifier.
        dstressdrateModifier = floatVector( currentStrain.size( ), 0 );

        //Compute the "infinite" term
        floatType scalarTerm = materialParameters[ 0 ];
        floatType taui, Gi, factor;
        std::vector< unsigned int > indices( currentStrain.size( ), 0 );
        floatVector Xic;
        unsigned int nTerms = ( materialParameters.size( ) - 1 )/2;
        unsigned int dim = currentStrain.size( );

        //Add the contributions from the other terms
        for ( unsigned int i=1; i<nTerms+1; i++ ){
            taui = materialParameters[ i ];
            Gi = materialParameters[ i+nTerms ];

            //Set the indices of the current values of the state variables
            for ( unsigned int j=dim*( i-1 ), k=0; j<dim*i; j++, k++ ){
                indices[ k ] = j;
            }

            //Get the previous values of the state variables
            tardigradeVectorTools::getValuesByIndex( currentStateVariables, indices, Xic );

            factor = 1./( taui + dt*( 1 - alpha )*currentRateModifier );

            //Add terms to the gradient w.r.t. the strain
            scalarTerm += Gi*( 1 - factor*( 1 - alpha )*dt*currentRateModifier );

            //Add terms to the gradient w.r.t. the rate modifier
            dstressdrateModifier -= Gi*( dt*( 1-alpha )/( taui + dt*( 1-alpha )*currentRateModifier ) )*( currentStrain - Xic );
        }

        //Assemble the full gradient
        dstressdstrain = scalarTerm*tardigradeVectorTools::eye<floatType>( currentStrain.size( ) );

        return;
    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables, floatMatrix &dstressdstrain,
                                   floatVector &dstressdrateModifier ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &previousStrain: The previous value of strain
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &stress: The computed stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         * \param &dstressdstrain: The derivative of the stress w.r.t. the strain
         * \param &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
         */

        floatVector dStress;

        return linearViscoelasticity( currentTime, currentStrain, previousTime, previousStrain,
                                     currentRateModifier, previousRateModifier,
                                     previousStateVariables, materialParameters, alpha,
                                     dStress, stress, currentStateVariables, dstressdstrain, dstressdrateModifier );

    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                   const floatType &previousTime, const floatVector &previousStrain,
                                   const floatType &currentRateModifier, const floatType &previousRateModifier,
                                   const floatVector &previousStateVariables, const floatVector &materialParameters,
                                   const floatType &alpha,
                                   floatVector &stress, floatVector &currentStateVariables,
                                   floatMatrix &dstressdstrain, floatVector &dstressdrateModifier,
                                   floatMatrix &dstressdPreviousStrain, floatVector &dstressdPreviousRateModifier,
                                   floatMatrix &dstressdPreviousStateVariables,
                                   floatMatrix &dStateVariablesdStrain, floatVector &dStateVariablesdRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStrain, floatVector &dStateVariablesdPreviousRateModifier,
                                   floatMatrix &dStateVariablesdPreviousStateVariables ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &previousStrain: The previous value of strain
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &stress: The computed stress in the reference configuration ( i.e. the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         * \param &dstressdstrain: The derivative of the stress w.r.t. the strain
         * \param &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
         * \param &dstressdPreviousStrain: The derivative of the stress w.r.t. the previous strain
         * \param &dstressdPreviousRateModifier: The derivative of the stress w.r.t. the previous rate modifier
         * \param &dstressdPreviousStateVariables: The derivative of the stress w.r.t. the previous state variables
         * \param &dStateVariablesdStrain: The derivative of the state variables w.r.t. the strain
         * \param &dStateVariablesdRateModifier: The derivative of the state variables w.r.t. the rate modifier
         * \param &dStateVariablesdPreviousStrain: The derivative of the state variables w.r.t. the previous strain
         * \param &dStateVariablesdPreviousRateModifier: The derivative of the state variables w.r.t. the previous rate modifier
         * \param &dStateVariablesdPreviousStateVariables: The derivative of the state variables w.r.t. the previous state variables
         */

        floatVector dStress;

        return linearViscoelasticity( currentTime, currentStrain, previousTime, previousStrain,
                                      currentRateModifier, previousRateModifier, previousStateVariables,
                                      materialParameters, alpha,
                                      dStress, stress, currentStateVariables,
                                      dstressdstrain, dstressdrateModifier, dstressdPreviousStrain, dstressdPreviousRateModifier,
                                      dstressdPreviousStateVariables,
                                      dStateVariablesdStrain, dStateVariablesdRateModifier,
                                      dStateVariablesdPreviousStrain, dStateVariablesdPreviousRateModifier,
                                      dStateVariablesdPreviousStateVariables );

    }

    void linearViscoelasticity( const floatType &currentTime, const floatVector &currentStrain,
                                    const floatType &previousTime, const floatVector &previousStrain,
                                    const floatType &currentRateModifier, const floatType &previousRateModifier,
                                    const floatVector &previousStateVariables, const floatVector &materialParameters,
                                    const floatType &alpha,
                                    floatVector &dStress, floatVector &stress, floatVector &currentStateVariables,
                                    floatMatrix &dstressdstrain, floatVector &dstressdrateModifier,
                                    floatMatrix &dstressdPreviousStrain, floatVector &dstressdPreviousRateModifier,
                                    floatMatrix &dstressdPreviousStateVariables,
                                    floatMatrix &dStateVariablesdStrain, floatVector &dStateVariablesdRateModifier,
                                    floatMatrix &dStateVariablesdPreviousStrain, floatVector &dStateVariablesdPreviousRateModifier,
                                    floatMatrix &dStateVariablesdPreviousStateVariables ){
        /*!
         * Compute the stress for linear viscoelasticity based on the potential function
         *
         * \f$\rho^0 \Psi = 0.5*( E_{ IJ } G_{ \infty } E_{ IJ } + \sum_{ n=1 }^N ( E_{ IJ } - \Xi_{ IJ }^n ) G^n ( E_{ IJ } - \Xi_{ IJ } ) )\f$
         *
         * \param &currentTime: The current time
         * \param &currentStrain: The current Green-Lagrange strain
         * \param &previousTime: The previous time
         * \param &previousStrain: The previous value of strain
         * \param &currentRateModifier: The current value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousRateModifier: The previous value of the rate modifier
         *     which can be used for temperature effects or ( potentially ) other non-linear effects.
         * \param &previousStateVariables: The previous values of the state variables
         * \param &materialParameters: The material parameters
         *     The order of the parameters is [ \f$G_{ \infty }\f$, \f$\tau\f$ s, \f$G\f$ s ] where
         *         - \f$G_{ \infty }\f$: The infinite stiffness modulus
         *         - \f$\tau\f$ s: The time constants
         *         - \f$G\f$ s: The stiffness values
         * \param &alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param &dStress: The change in stress in the reference configuration ( i.e., the same configuration as the strain )
         * \param &stress: The computed stress in the reference configuration ( i.e., the same configuration as the strain )
         * \param &currentStateVariables: The current values of the state variables
         * \param &dstressdstrain: The derivative of the stress w.r.t. the strain
         * \param &dstressdrateModifier: The derivative of the stress w.r.t. the rate modifier
         * \param &dstressdPreviousStrain: The derivative of the stress w.r.t. the previous strain
         * \param &dstressdPreviousRateModifier: The derivative of the stress w.r.t. the previous rate modifier
         * \param &dstressdPreviousStateVariables: The derivative of the stress w.r.t. the previous state variables
         * \param &dStateVariablesdStrain: The derivative of the state variables w.r.t. the strain
         * \param &dStateVariablesdRateModifier: The derivative of the state variables w.r.t. the rate modifier
         * \param &dStateVariablesdPreviousStrain: The derivative of the state variables w.r.t. the previous strain
         * \param &dStateVariablesdPreviousRateModifier: The derivative of the state variables w.r.t. the previous rate modifier
         * \param &dStateVariablesdPreviousStateVariables: The derivative of the state variables w.r.t. the previous state variables
         */
        floatType dt = currentTime - previousTime;

        //Check the material parameters
        TARDIGRADE_ERROR_TOOLS_CHECK( materialParameters.size( ) != 0, "At least one material parameter must be provided." )

        //Set the number of Prony-Series terms
        TARDIGRADE_ERROR_TOOLS_CHECK( ( materialParameters.size( ) - 1 ) % 2 == 0, "An equal number of taus and Gs must be provided." )

        unsigned int nTerms = ( materialParameters.size( ) - 1 ) / 2;

        //Set the dimension of the strain
        unsigned int dim = currentStrain.size( );
        //Check the state variables
        TARDIGRADE_ERROR_TOOLS_CHECK( previousStateVariables.size( ) == nTerms * dim, "The number of previous state variables is not consistent with the strain size." )

        //Compute the infinite stress
        floatVector dStrain = currentStrain - previousStrain;
        dStress = materialParameters[ 0 ]*dStrain;
        stress = materialParameters[ 0 ]*currentStrain;

        floatMatrix EYE = tardigradeVectorTools::eye< floatType >( dim );

        dstressdstrain = materialParameters[ 0 ] * EYE;
        dstressdPreviousStrain = floatMatrix( stress.size( ), floatVector( previousStrain.size( ), 0 ) );
        dstressdPreviousStateVariables = floatMatrix( stress.size( ), floatVector( previousStateVariables.size( ), 0 ) );

        //Set the initial value of the factor
        floatType factor;
        floatType dfactordr;
        floatType taui;
        floatType Gi;
        floatVector Xip( dim, 0 ), Xic( dim, 0 );
        currentStateVariables.resize( 0 );

        std::vector< unsigned int > indices( dim, 0 );

        floatVector dXi( dim, 0 );

        dstressdrateModifier = floatVector( currentStrain.size( ), 0 );
        dstressdPreviousRateModifier = floatVector( currentStrain.size( ), 0 );

        dStateVariablesdStrain = floatMatrix( previousStateVariables.size( ), floatVector( currentStrain.size( ), 0 ) );
        dStateVariablesdPreviousStrain = floatMatrix( previousStateVariables.size( ), floatVector( previousStrain.size( ), 0 ) );

        floatMatrix dStateVariablesdCurrentRateModifier_matrix( nTerms );
        floatMatrix dStateVariablesdPreviousRateModifier_matrix( nTerms );
        dStateVariablesdPreviousStateVariables = floatMatrix( previousStateVariables.size( ), floatVector( previousStateVariables.size( ), 0 ) );

        floatType dstressdstrain_relax = 0;
        floatType dstressdPreviousStrain_relax = 0;

        for ( unsigned int i = 1; i < nTerms + 1; i++ ){
            //Get the parameters
            taui = materialParameters[ i ];
            Gi = materialParameters[ i + nTerms ];

            //Get the factor
            factor = dt/( taui + dt*( 1 - alpha )*currentRateModifier );
            dfactordr = -( dt * dt ) * ( 1 - alpha ) / std::pow( ( taui + dt * ( 1 - alpha ) * currentRateModifier ), 2 );

            //Set the indices of the previous values of the state variables
            for ( unsigned int j = dim * ( i - 1 ), k = 0; j < dim * i; j++, k++ ){

                indices[ k ] = j;

            }

            //Get the previous values of the state variables
            tardigradeVectorTools::getValuesByIndex( previousStateVariables, indices, Xip );

            //Compute the new state-variable values
            floatVector dxi = factor * ( ( 1 - alpha ) * ( currentStrain - Xip ) * currentRateModifier
                                             + alpha * ( previousStrain - Xip ) * previousRateModifier );

            dStateVariablesdCurrentRateModifier_matrix[ i - 1 ] = dfactordr * ( ( 1 - alpha ) * ( currentStrain - Xip ) * currentRateModifier
                                                                                  + alpha * ( previousStrain - Xip ) * previousRateModifier )
                                                                + factor * ( 1 - alpha ) * ( currentStrain - Xip );

            dStateVariablesdPreviousRateModifier_matrix[ i - 1 ] = factor * alpha * ( previousStrain - Xip );

            floatType dxidcurrentStrain_factor = factor * ( 1 - alpha ) * currentRateModifier;

            floatType dxidpreviousStrain_factor = factor * alpha * previousRateModifier;

            floatType dxidPreviousXi_factor = -factor * ( ( 1 - alpha ) * currentRateModifier
                                                              + alpha   * previousRateModifier );

            Xic = Xip + dxi;

            //Add the contribution to the stress
            dStress += Gi * ( dStrain - dxi );

            stress += Gi * ( currentStrain - Xic );

            dstressdstrain_relax += Gi * ( 1 - dxidcurrentStrain_factor );

            dstressdrateModifier += -Gi * dStateVariablesdCurrentRateModifier_matrix[ i - 1 ];

            dstressdPreviousStrain_relax += -Gi * dxidpreviousStrain_factor;

            dstressdPreviousRateModifier += -Gi * dStateVariablesdPreviousRateModifier_matrix[ i - 1 ];

            for ( unsigned int j = 0; j < dim; j++ ){

                dstressdPreviousStateVariables[ j ][ dim * ( i - 1 ) + j ] = -Gi * ( 1 + dxidPreviousXi_factor );

                dStateVariablesdStrain[ dim * ( i - 1 ) + j ][ j ] = dxidcurrentStrain_factor;

                dStateVariablesdPreviousStrain[ dim * ( i - 1 ) + j ][ j ] = dxidpreviousStrain_factor;

                dStateVariablesdPreviousStateVariables[ dim * ( i - 1 ) + j ][ dim * ( i - 1 ) + j ]
                    = 1 + dxidPreviousXi_factor;

            }

            //Save the new value of the state variable
            currentStateVariables = tardigradeVectorTools::appendVectors( { currentStateVariables, Xic } );
        }

        dstressdstrain += dstressdstrain_relax * EYE;

        dstressdPreviousStrain += dstressdPreviousStrain_relax * EYE;

        dStateVariablesdRateModifier = tardigradeVectorTools::appendVectors( dStateVariablesdCurrentRateModifier_matrix );

        dStateVariablesdPreviousRateModifier = tardigradeVectorTools::appendVectors( dStateVariablesdPreviousRateModifier_matrix );

        return;

    }

    void volumetricNeoHookean( const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         *
         * \f$U( J ) = 0.5*bulkModulus*( 0.5*( J**2 - 1 ) - ln( J ) )\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &jacobian: The jacobian of deformation
         * \param &bulkModulus: The bulk modulus
         * \param &meanStress: The mean stress of the material. NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         *
         */

        //Error handling
        TARDIGRADE_ERROR_TOOLS_CHECK( jacobian>0, "determinant is less than or equal zero" );

        //Compute the meanStress
        meanStress = 0.5*bulkModulus*( jacobian - 1/jacobian );

        return;
    }

    void volumetricNeoHookean( const floatType &jacobian, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         *
         * \f$U( J ) = 0.5*bulkModulus*( 0.5*( J**2 - 1 ) - ln( J ) )\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &jacobian: The jacobian of deformation
         * \param &bulkModulus: The bulk modulus
         * \param &meanStress: The mean stress of the material. NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         * \param &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian
         *     of deformation.
         */

        //Error handling
        TARDIGRADE_ERROR_TOOLS_CHECK( jacobian>0, "determinant is less than or equal zero" )

        //Compute the meanStress
        meanStress = 0.5*bulkModulus*( jacobian - 1/jacobian );

        //Compute the derivative of the meanStress w.r.t. jacobian
        dmeanStressdJ = 0.5*bulkModulus*( 1 + 1/( jacobian*jacobian ) );

        return;
    }


    void volumetricNeoHookean( const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         *
         * \f$U( J ) = 0.5*bulkModulus*( 0.5*( J**2 - 1 ) - ln( J ) )\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &bulkModulus: The bulk modulus
         * \param &meanStress: The mean stress of the material. NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         *
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        //Check the size of the deformation gradient
        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size( ) == sot_dim, "deformation gradient must have nine terms" )

        //Compute the determinant of the deformation gradient
        floatType J;
        J = tardigradeVectorTools::determinant( deformationGradient, dim, dim );

        return volumetricNeoHookean( J, bulkModulus, meanStress );

    }

    void volumetricNeoHookean( const floatVector &deformationGradient, const floatType &bulkModulus,
                                  floatType &meanStress, floatType &dmeanStressdJ ){
        /*!
         * Compute the volumetric part of a Neo-Hookean material model response of the form
         *
         * \f$U( J ) = 0.5*bulkModulus*( 0.5*( J**2 - 1 ) - ln( J ) )\f$
         *
         * where \f$J\f$ is the determinant of the deformation gradient.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &bulkModulus: The bulk modulus
         * \param &meanStress: The mean stress of the material. NOTE: It is up to the user to determine
         *     which configuration this mean stress is defined within. If it is the reference
         *     configuration a mapping may be necessary going into a different configuration.
         * \param &dmeanStressdJ: The derivative of the mean stress w.r.t. the jacobian
         *     of deformation.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        //Check the size of the deformation gradient
        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size( ) == sot_dim, "deformation gradient must have nine terms" )

        //Compute the determinant of the deformation gradient
        floatType J;
        J = tardigradeVectorTools::determinant( deformationGradient, 3, 3 );

        return volumetricNeoHookean( J, bulkModulus, meanStress, dmeanStressdJ );

    }

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatType n, floatType &p ){
        /*!
         * Implementation of the Perzyna type model of the form
         *
         * \f$p = A \left \langle \frac{ f }{ q } \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param n: The exponent.
         * \param &p: The value of the model.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( !tardigradeVectorTools::fuzzyEquals( q, 0. ), "The denominator term is zero" )

        TARDIGRADE_ERROR_TOOLS_CHECK( n >= 1, "n must be >= 1" )

        p = A*pow( tardigradeConstitutiveTools::mac( f/q ), n );
        return;
    }

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p ){
        /*!
         * Implementation of the Perzyna type model of the form
         *
         * \f$p = A \left \langle \frac{ f }{ q } \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param &parameters: The parameters ( currently just n )
         * \param &p: The value of the model.
         */
        TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == 1, "The parameters vector is one value long" )

        return perzynaModel( f, q, A, parameters[ 0 ], p );
    }

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatType n, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA ){

        /*!
         * Implementation of the Perzyna type model of the form
         *
         * \f$p = A \left \langle \frac{ f }{ q } \right \rangle^n\f$
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

        TARDIGRADE_ERROR_TOOLS_CHECK( !tardigradeVectorTools::fuzzyEquals( q, 0. ), "The denominator term is zero" )

        TARDIGRADE_ERROR_TOOLS_CHECK( n >= 1, "n must be >= 1" )

        //Compute the value
        floatType mac, dmacdx;
        mac = tardigradeConstitutiveTools::mac( f/q, dmacdx );

        p = A*pow( tardigradeConstitutiveTools::mac( f/q ), n );

        dpdf = A*n*pow( mac, n-1 )*dmacdx/q;
        dpdq = -A*n*pow( mac, n-1 )*dmacdx*f/( q*q );
        dpdA = pow( tardigradeConstitutiveTools::mac( f/q ), n );

        return;
    }

    void perzynaModel( const floatType f, const floatType q, const floatType A, const floatVector &parameters, floatType &p,
                          floatType &dpdf, floatType &dpdq, floatType &dpdA ){
        /*!
         * Implementation of the Perzyna type model of the form
         *
         * \f$p = A \left \langle \frac{ f }{ q } \right \rangle^n\f$
         *
         * where \f$\langle\f$ \f$\rangle\f$ are the Macaulay brackets.
         *
         * \param f: The numerator term in the brackets.
         * \param q: The denominator term in the brackets.
         * \param A: The scaling factor.
         * \param &parameters: The parameters ( currently just n )
         * \param &p: The value of the model.
         * \param &dpdf: The derivative of the value w.r.t. f.
         * \param &dpdq: The derivative of the value w.r.t. q.
         * \param &dpdA: The derivative of the value w.r.t. A.
         */
        TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == 1, "The parameters vector is one value long" )

        return perzynaModel( f, q, A, parameters[ 0 ], p, dpdf, dpdq, dpdA );
    }

    void linearHardening( const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value ){
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

        TARDIGRADE_ERROR_TOOLS_CHECK( stateVariables.size( ) == linearModuli.size( ), "The state variables and the moduli must have the same size" );

        value = tardigradeVectorTools::dot( stateVariables, linearModuli ) + scalarShift;
        return;
    }

    void linearHardening( const floatVector &stateVariables, const floatVector &linearModuli, const floatType &scalarShift,
                             floatType &value, floatVector &valueJacobian ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( linearHardening( stateVariables, linearModuli, scalarShift, value ) )

        valueJacobian = linearModuli;
        return;
    }

    void computeJaumannStiffnessTensor( const floatVector &cauchyStress, const floatVector &currentDeformationGradient,
                                            const floatMatrix &dCauchydF, floatMatrix &C ){
        /*!
         * Compute the Jaumann stiffness tensor from the cauchy stress, the current deformation gradient,
         * and the total derivative of the Cauchy stress w.r.t. the deformation gradient
         * 
         * From the variation of the Jaumann rate
         * \f$J \mathbb{ C }_{ ijkl } \delta D_{ kl } = \delta \left( J \sigma_{ ij }\right) + \delta W_{ ik } \sigma_{ kj } - \sigma_{ ik } \delta W_{ kj }\f$
         *
         * Where \f$J$\f$ is the Jacobian of deformation, \f$\mathbb{ C }\f$ is the Jaumann stiffness tensor,
         * \f$\bf{ \sigma }\f$ is the Cauchy stress and \f$\delta \bf{ D }\f$ is the perturbation of the rate of
         * deformation, and \f$\delta \bf{ W }\f$ is the perturbation of the rate of spin.
         *
         * By using the properties that
         *
         * \f$\delta D_{ kl } = \text{ symm }\left( \delta F_{ kK } F_{ Kl }^{ -1 }\right)\f$
         *
         * \f$\delta W_{ kl } = \text{ asymm }\left( \delta F_{ kK } F_{ Kl }^{ -1 }\right)\f$
         *
         * Where
         *
         * \f$\text{ symm }\left( \bf{ A }\right ) = \frac{ 1 }{ 2 }\left( \bf{ A } + \bf{ A }^T\right)\f$
         *
         * \f$\text{ asymm }\left( \bf{ A }\right ) = \frac{ 1 }{ 2 }\left( \bf{ A } - \bf{ A }^T\right)\f$
         *
         * It can be shown that
         *
         * \f$\mathbb{ C }_{ ijkl } = \left( \delta_{ mn } \sigma_{ ij } + \frac{ D \sigma_{ ij } }{ D F_{ mK } } F_{ nK }\right)\mathbb{ P }_{ mnkl }^{ symm }\f$
         *
         * Where
         *
         * \f$\mathbb{ P }_{ ijkl }^{ symm } = \frac{ 1 }{ 2 }\left( \delta_{ ik } \delta_{ jl } + \delta_{ jk } \delta_{ il }\right)\f$
         * 
         * \param &cauchyStress: The Cauchy stress
         * \param &currentDeformationGradient: The current deformation gradient i.e. the mapping for differential
         *     lengths from the reference to the current configuration
         * \param &dCauchydF: The Jacobian ( total derivative ) of the Cauchy stress w.r.t. the deformation
         *     gradient
         * \param &C: The Jaumann stiffness tensor
         */

        // Perform error handling
        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == currentDeformationGradient.size( ), "The cauchy stress ( length " + std::to_string( cauchyStress.size( ) ) + " ) and the deformation gradient ( length " + std::to_string( currentDeformationGradient.size( ) ) + " ) must have the same size" )

        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == dCauchydF.size( ), "The derivative of the Cauchy stress w.r.t. the deformation gradient has " + std::to_string( dCauchydF.size( ) ) + " rows. It needs to have " + std::to_string( cauchyStress.size( ) ) + " to be consistent with the provided Cauchy stress" );

        for ( unsigned int i = 0; i < dCauchydF.size( ); i++ ){

            TARDIGRADE_ERROR_TOOLS_CHECK( dCauchydF[ i ].size( ) == currentDeformationGradient.size( ), "Row " + std::to_string( i ) + " of dCauchydF is of length " + std::to_string( dCauchydF[ i ].size( ) ) + " and it should have a length of " + std::to_string( currentDeformationGradient.size( ) ) )

        }

        // Build the second order identity tensor
        floatVector eye( cauchyStress.size( ), 0 );
        tardigradeVectorTools::eye( eye );

        // Set the dimension of the problem
        unsigned int dim = tardigradeVectorTools::trace( eye );

        // Construct the symmetric projection tensor
        floatMatrix Psymm( cauchyStress.size( ), floatVector( cauchyStress.size( ), 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    for ( unsigned int l = 0; l < dim; l++ ){

                        Psymm[ dim * i + j ][ dim * k + l ] =
                            0.5 * ( eye[ dim * i + k ] * eye[ dim * j + l ] + eye[ dim * j + k ] * eye[ dim * i + l ] );

                    }

                }

            }

        }

        // Initialize the stiffness tensor by including cauchyStress dyad eye
        C = tardigradeVectorTools::dyadic( cauchyStress, eye );

        // Add the dCauchydF term

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int m = 0; m < dim; m++ ){

                    for ( unsigned int n = 0; n < dim; n++ ){

                        for ( unsigned int K = 0; K < dim; K++ ){

                                C[ dim * i + j ][ dim * m + n ] +=
                                    dCauchydF[ dim * i + j ][ dim * m + K ] * currentDeformationGradient[ dim * n + K ];

                        }

                    }

                }

            }

        }

        // Get the symmetric part
        C = tardigradeVectorTools::dot( C, Psymm );

        return;

    }

}
