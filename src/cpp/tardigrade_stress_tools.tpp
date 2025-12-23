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
    void TARDIGRADE_OPTIONAL_INLINE calculateMeanStress(const stress_iterator &stress_begin,
                                                        const stress_iterator &stress_end, stress_type &meanStress);

    template <unsigned int dim, typename stress_type, class stress_iterator>
    stress_type TARDIGRADE_OPTIONAL_INLINE calculateMeanStress(const stress_iterator &stress_begin,
                                                               const stress_iterator &stress_end);

    template <unsigned int dim, class stress_iterator, typename stress_type, class jacobian_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateMeanStress(const stress_iterator &stress_begin,
                                                        const stress_iterator &stress_end, stress_type &meanStress,
                                                        jacobian_iterator jacobian_begin,
                                                        jacobian_iterator jacobian_end);

    template <unsigned int dim, class stress_iterator, class deviatoric_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(const stress_iterator &stress_begin,
                                                              const stress_iterator &stress_end,
                                                              deviatoric_iterator    deviatoric_begin,
                                                              deviatoric_iterator    deviatoric_end);

    template <unsigned int dim, class stress_iterator, class deviatoric_iterator, class jacobian_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateDeviatoricStress(
        const stress_iterator &stress_begin, const stress_iterator &stress_end, deviatoric_iterator deviatoric_begin,
        deviatoric_iterator deviatoric_end, jacobian_iterator jacobian_begin, jacobian_iterator jacobian_end);

    template <unsigned int dim, class stress_iterator, typename vonMises_type>
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(const stress_iterator &stress_begin,
                                                            const stress_iterator &stress_end, vonMises_type &vonMises);

    template <unsigned int dim, class stress_iterator, typename vonMises_type, class jacobian_iterator>
    void TARDIGRADE_OPTIONAL_INLINE calculateVonMisesStress(const stress_iterator &stress_begin,
                                                            const stress_iterator &stress_end, vonMises_type &vonMises,
                                                            jacobian_iterator jacobian_begin,
                                                            jacobian_iterator jacobian_end);

}  // namespace tardigradeStressTools
