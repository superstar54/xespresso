.. _nscf:

===========================================
nscf
===========================================


A example of nscf calculation following the above one.

.. code-block:: python

    from xespresso.post.nscf import EspressoNscf
    nscf = EspressoNscf(calc.directory, prefix = calc.prefix,
                    occupations = 'tetrahedra',
                    kpts = (2, 2, 2),
                    debug = True,
                    )
    nscf.run()
