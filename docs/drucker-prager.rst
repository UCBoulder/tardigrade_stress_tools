Drucker-Prager Yield Criterion
==============================

The `Drucker-Prager yield criterion
<https://en.wikipedia.org/wiki/Drucker%E2%80%93Prager_yield_criterion>`_ is
defined as 

.. math::

    0 = \sigma^{vonMises} - A \sigma^{mean} - B
    :label: drucker-prager

For use in visco-plastic flow equations, Equation :eq:`drucker-prager` may be
re-written as

.. math::

    f = \sigma^{vonMises} - A \sigma^{mean} - B
    :label: drucker-prager-flow

Equation :eq:`drucker-prager-flow` is implemented as ``druckerPragerSurface``.
In strain-hardening materials and visco-plastic materials it is possible to
cross the initial yield surface, where ``f >= 0`` indicates inelastic flow, ``f
= 0`` is inelastic flow on the 'yield surface', and ``f <= 0`` indicates no (or
minimal) inelastic flow. 

For these cases, it is important to also compute the inelastic strain increment
direction. For associative flow this may be written as

.. math::

    n = \frac{\partial f}{\partial \sigma}  / \norm{\frac{\partial f}{\partial
        \sigma}}
    :label: drucker-prager-direction

and for non-associative flow

.. math::

    g = \sigma^{vonMises} - C \sigma^{mean} - D

    n = \frac{\partial g}{\partial \sigma}  / \norm{\frac{\partial g}{\partial
        \sigma}}

where

.. math::

    \frac{\partial g}{\partial \sigma} = \frac{\partial
        \sigma^{vonMises}}{\sigma} - C \frac{\sigma^{mean}}{\sigma}
    :label: drucker-prager-non-associative

Equation :eq:`drucker-prager-non-associative` is implemented as
``druckerPragerDirection`` and the derivation of the R.H.S. partial derivatives
are implemented in ``vonMises`` and ``meanStress``, respectively. 

The von Mises stress is defined as

.. math::

    \sigma^{vonMises} = \sqrt{\frac{3}{2} S S}
    \sigma^{vonMises} = \sqrt{\frac{3}{2} S_{ij} S_{ij}}

where ``S`` is the deviatoric stress

.. math::

    S = \sigma - \frac{1}{3} \sigma \delta
    S_{ij} = \sigma_{ij} - \frac{1}{3} \sigma_{ij} \delta_{ij}

and ``\delta`` is the `Kronecker delta
<https://en.wikipedia.org/wiki/Kronecker_delta>`_. The mean stress is defined as

.. math::

    \sigma^{mean} = \frac{1}{3} trace \left ( \sigma \right )
    \sigma^{mean} = \frac{1}{3} trace \left ( \sigma_{ii} \right )
 
From these definitions, we can calculate the partial derivatives in Equation
:eq:`drucker-prager-non-associative`. First, re-write Equation
:eq:`drucker-prager-non-associative` with the chain rule.

.. math::

    \frac{\partial g}{\partial \sigma_{kl}} = \frac{\partial
        \sigma^{vonMises}}{\S_{ij}}\frac{\partial \S_{ij}}{\partial \sigma_{kl}}
        - C \frac{\sigma^{mean}}{\sigma_{kl}} 
        - \frac{\partial D}{\partial \sigma_{kl}}
    :label: drucker-prager-chain-rule
    
The partial derivative of the constant ``D`` is zero. The partial derivative of
von Mises stres with respect to the deviatoric stress tensor can be computed
with the property of `second-order tensor derivatives
<https://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)>`_

.. math::

    \frac{\partial x_{ij}}{\partial y_{kl}} = \delta_{ik} \delta_{jl}
    :label: tensor-partial

and the property of the Kronecker delta

.. math::

    x_{kl} = x_{ij} \delta_{ik} \delta_{jl}
    :label: kronecker-property

First, compute the partial derivative of the von Mises stress with respect to
the deviatoric stress.

.. math::

    \frac{\partial \sigma^{vonMises}}{S_{kl}} = \frac{1}{2 \sigma^{vonMises}}
        \frac{3 S_{ij}}{2} \left ( \delta_{ik} \delta_{jl} + \delta_{ik}
        \delta_{jl}

    \frac{\partial \sigma^{vonMises}}{S_{kl}} = \frac{1}{2 \sigma^{vonMises}
        \frac{3 S_{ij}}{2} 2 \delta_{ik} \delta_{jl}

    \frac{\partial \sigma^{vonMises}}{S_{kl}} = \frac{3}{2 \sigma^{vonMises}}
        \S_{ij} \delta_{ik} \delta_{jl}

    \frac{\partial \sigma^{vonMises}}{S_{kl}} = \frac{3 S_{ij}}{2
        \sigma^{vonMises}}
    
Next, compute the partial derivative of the mean stress with respect to the
stress tensor. 

.. math::

    \frac{\partial \sigma^{mean}}{\partial \sigma_{kl}} = \frac{1}{3}
        \delta_{ik} \delta_{il}

    \frac{\partial \sigma^{mean}}{\partial \sigma_{kl}} = \frac{1}{3} \delta_{kl}

Finally, compute partial derivative of the deviatoric stress with respect to the
stress. 

.. math::

    \frac{\partial S_{ij}}{\partial \sigma_{kl}} = \delta_{ik} \delta_{jl} -
        \frac{\partial \sigma^{mean}}{\partial \sigma_{kl}} \delta_{ij}

    \frac{\partial S_{ij}}{\partial \sigma_{kl}} = \delta_{ik} \delta_{jl} -
        \frac{1}{3} \delta_{kl} \delta_{ij}

These partial derivatives are implemented in ``calculateVonMises``,
``calculateDeviatoricStress``, and ``calculateMeanStress`` and are used in
``druckerPragerSurface`` to calculate the flow direction.

For completeness, the Drucker-Prager flow direction is included below.

.. math::

    
    \frac{\partial g}{\partial \sigma_{kl}} = \frac{\partial
        \sigma^{vonMises}}{\S_{ij}}\frac{\partial \S_{ij}}{\partial \sigma_{kl}}
        - C \frac{\sigma^{mean}}{\sigma_{kl}} 
        - \frac{\partial D}{\partial \sigma_{kl}}

    \frac{\partial g}{\partial \sigma_{kl}} = \frac{3 S_{ij}}{2
        \sigma^{vonMises}} \left ( \delta_{ik} \delta_{jl} -
        \frac{1}{3} \delta_{kl} \delta_{ij} \right ) - C \frac{1}{3} \delta_{kl}

The solution may be simplified because multiplication of the deviatoric stress
by the Kronecker delta with matching dimensions results in the zero valued
tensor. 

.. math::

    0 = S_{ij} \delta_{ij}

Simplifying, the final result for the partial derivative of the flow direction
is

.. math::
 
    \frac{\partial g}{\partial \sigma_{kl}} = \frac{3 S_{kl}}{2
        \sigma^{vonMises}} - \frac{A}{3} \delta_{kl}
