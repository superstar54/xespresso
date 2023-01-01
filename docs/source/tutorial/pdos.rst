.. _pdos:

===========================================
pdos
===========================================


A example of calculating and plotting the pdos from the nscf calculation.

.. code-block:: python

    from xespresso.post.dos import EspressoDos
    # pdos
    from xespresso.post.projwfc import EspressoProjwfc
    projwfc = EspressoProjwfc(parent_directory = 'calculations/scf/co',
                prefix = 'co',
                DeltaE = 0.01)
    projwfc.run()

.. image:: /_static/images/co-pdos.png
   :width: 10cm

List of all Methods
--------------------

.. autoclass:: xespresso.post.projwfc.EspressoProjwfc
   :members:
