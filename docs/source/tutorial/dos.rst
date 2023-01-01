.. _dos:

===========================================
dos
===========================================


A example of calculating and plotting the pdos from the nscf calculation.

.. code-block:: python

    from xespresso.post.dos import EspressoDos
    # dos
    dos = EspressoDos(parent_directory = 'calculations/scf/co',
                prefix = calc.prefix,
                Emin = fe - 30, Emax = fe + 30, DeltaE = 0.01)
    dos.run()


.. image:: /_static/images/co-dos.png
   :width: 10cm


List of all Methods
--------------------

.. autoclass:: xespresso.post.dos.EspressoDos
   :members:
