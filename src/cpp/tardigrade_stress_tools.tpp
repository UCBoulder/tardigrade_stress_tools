/**
 ******************************************************************************
 * \file tardigrade_stress_tools.tpp
 ******************************************************************************
 * A collection of tools which implement and solve stress-strain relationships
 * in such a way to enable more rapid development of constitutive models which
 * have capabilities which may not be contained within a collection of
 * of constitutive models.
 ******************************************************************************
 */

namespace tardigradeStressTools {

    template <unsigned int dim, class stress_iterator, typename stress_type>
    void calculateMeanStress(const stress_iterator &stress_begin, const stress_iterator &stress_end,
                             stress_type &meanStress) {
        /*!
         * Compute the mean stress from a 2nd rank stress tensor stored in row major format
         *
         * \f$meanStress = \frac{ 1 }{ 3 }*trace( \sigma )\f$
         *
         * \param &stress_begin: The starting iterator of the stress tensor
         * \param &stress_end: The stopping iterator of the stress tensor
         * \param &meanStress: The mean stress scalar
         */

        tardigradeVectorTools::rowMajorTrace<dim, dim>(stress_begin, stress_end, meanStress);

        meanStress /= 3.;
    }

    template <unsigned int dim, typename stress_type, class stress_iterator>
    stress_type calculateMeanStress(const stress_iterator &stress_begin, const stress_iterator &stress_end) {
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
        calculateMeanStress<dim>(stress_begin, stress_end, meanStress);

        return meanStress;
    }

    template <unsigned int dim, class stress_iterator, typename stress_type, class jacobian_iterator>
    void calculateMeanStress(const stress_iterator &stress_begin, const stress_iterator &stress_end,
                             stress_type &meanStress, jacobian_iterator jacobian_begin,
                             jacobian_iterator jacobian_end) {
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

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(jacobian_end - jacobian_begin) == dim * dim,
                                     "jacobian has a size of " + std::to_string(jacobian_end - jacobian_begin) +
                                         " but should have a size of " + std::to_string(dim * dim));

        calculateMeanStress<dim>(stress_begin, stress_end, meanStress);

        std::fill(jacobian_begin, jacobian_end, jacobian_type());

        for (unsigned int i = 0; i < dim; ++i) {
            *(jacobian_begin + dim * i + i) = 1. / 3.;
        }
    }

    template <unsigned int dim, class stress_iterator, class deviatoric_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(const stress_iterator &stress_begin,
                                                              const stress_iterator &stress_end,
                                                              deviatoric_iterator    deviatoric_begin,
                                                              deviatoric_iterator    deviatoric_end) {
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

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(stress_end - stress_begin) == dim * dim,
                                     "The stress has a size of " +
                                         std::to_string((unsigned int)(stress_end - stress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deviatoric_end - deviatoric_begin) == dim * dim,
                                     "The deviatoric stress has a size of " +
                                         std::to_string((unsigned int)(deviatoric_end - deviatoric_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim));

        stress_type meanStress;

        calculateMeanStress<dim>(stress_begin, stress_end, meanStress);

        std::copy(stress_begin, stress_end, deviatoric_begin);

        for (unsigned int i = 0; i < dim; ++i) {
            *(deviatoric_begin + dim * i + i) -= meanStress;
        }

        return;
    }

    template <unsigned int dim, class stress_iterator, class deviatoric_iterator, class jacobian_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end, deviatoric_iterator deviatoric_begin,
        deviatoric_iterator deviatoric_end, jacobian_iterator jacobian_begin, jacobian_iterator jacobian_end) {
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

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(jacobian_end - jacobian_begin) == dim * dim * dim * dim,
                                     "The jacobian has a size of " +
                                         std::to_string((unsigned int)(jacobian_end - jacobian_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim * dim * dim));

        calculateDeviatoricStress<dim>(stress_begin, stress_end, deviatoric_begin, deviatoric_end);

        std::fill(jacobian_begin, jacobian_end, stress_type());

        for (unsigned int i = 0; i < dim * dim; ++i) {
            *(jacobian_begin + dim * dim * i + i) += 1;
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                *(jacobian_begin + dim * dim * dim * i + dim * dim * i + dim * j + j) -= 1. / 3;
            }
        }

        return;
    }

    template <unsigned int dim, class stress_iterator, typename vonMises_type>
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(const stress_iterator &stress_begin,
                                                            const stress_iterator &stress_end,
                                                            vonMises_type         &vonMises) {
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

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(stress_end - stress_begin) == dim * dim,
                                     "The stress tensor has a size of " +
                                         std::to_string((unsigned int)(stress_end - stress_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim));

        std::array<stress_type, dim * dim> deviatoric = {};
        calculateDeviatoricStress<dim>(stress_begin, stress_end, std::begin(deviatoric), std::end(deviatoric));

        vonMises = std::sqrt(1.5 * std::inner_product(std::begin(deviatoric), std::end(deviatoric),
                                                      std::begin(deviatoric), stress_type()));
    }

    template <unsigned int dim, class stress_iterator, typename vonMises_type, class jacobian_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(const stress_iterator &stress_begin,
                                                            const stress_iterator &stress_end, vonMises_type &vonMises,
                                                            jacobian_iterator jacobian_begin,
                                                            jacobian_iterator jacobian_end) {
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

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(stress_end - stress_begin) == dim * dim,
                                     "The stress tensor has a size of " +
                                         std::to_string((unsigned int)(stress_end - stress_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(stress_end - stress_begin) == dim * dim,
                                     "The jacobian has a size of " +
                                         std::to_string((unsigned int)(jacobian_end - jacobian_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim));

        std::array<stress_type, dim * dim> deviatoric = {};
        calculateDeviatoricStress<dim>(stress_begin, stress_end, std::begin(deviatoric), std::end(deviatoric));

        vonMises = std::sqrt(1.5 * std::inner_product(std::begin(deviatoric), std::end(deviatoric),
                                                      std::begin(deviatoric), stress_type()));

        tardigradeConstitutiveTools::computeUnitNormal(std::begin(deviatoric), std::end(deviatoric), jacobian_begin,
                                                       jacobian_end);

        std::transform(jacobian_begin, jacobian_end, jacobian_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, std::sqrt(1.5)));
    }

    template <typename vonMises_type, typename meanStress_type, class dpParam_iterator, typename dpYield_type>
    void druckerPragerSurface_iter(const vonMises_type &vonMises, const meanStress_type &meanStress,
                                   const dpParam_iterator &dpParam_begin, const dpParam_iterator &dpParam_end,
                                   dpYield_type &dpYield) {
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
         * \param &dpParam_begin: The starting iterator of the Drucker-Prager material parameters
         * \param &dpParam_end: The stopping iterator of the Drucker-Prager material parameters
         * \param &dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dpParam_end - dpParam_begin) == 2,
                                     "The Drucker-Prager parameter vector has a size of " +
                                         std::to_string((unsigned int)(dpParam_end - dpParam_begin)) +
                                         " but must have a size of 2");

        druckerPragerSurface(vonMises, meanStress, *(dpParam_begin + 0), *(dpParam_begin + 1), dpYield);
    }

    template <typename vonMises_type, typename meanStress_type, class dpParam_iterator, typename dpYield_type>
    dpYield_type druckerPragerSurface_iter(const vonMises_type &vonMises, const meanStress_type &meanStress,
                                           const dpParam_iterator &dpParam_begin, const dpParam_iterator &dpParam_end) {
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
         * \param &dpParam_begin: The starting iterator of the Drucker-Prager material parameters
         * \param &dpParam_end: The stopping iterator of the Drucker-Prager material parameters
         * \return dpYield: The Drucker-Prager yield stress/criterion/surface
         */

        dpYield_type dpYield;

        druckerPragerSurface_iter(vonMises, meanStress, dpParam_begin, dpParam_end, dpYield);

        return dpYield;
    }

}  // namespace tardigradeStressTools
