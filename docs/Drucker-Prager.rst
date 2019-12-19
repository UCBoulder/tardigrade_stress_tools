The Drucker-Prager yield criterion is defined as 

.. math::

    0 = \sigma^{vonMises} - A * \sigma^{mean} - B
    :label: drucker-prager

For use in visco-plastic flow equations, Equation :eq:`drucker-prager` may be
re-written as

.. math::

    f = \sigma^{vonMises} - A * \sigma^{mean} - B
    :label: drucker-prager-flow

Equation :eq:`drucker-prager-flow` is implemented as ``druckerPragerSurface``.
In strain-hardening materials and visco-plastic materials it is possible to
cross the initial yield surface. For these cases, it is important to also
compute the inelastic strain increment direction. For associative flow this may
be written as

.. math::

    n = \frac{\partial f}{\partial \sigma}

and for non-associative flow

.. math::

    f = \sigma^{vonMises} - C * \sigma^{mean} - D
    n = \frac{\partial g}{\partial \sigma} 
