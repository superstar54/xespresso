.. _neb:

===========================================
NEB
===========================================
A example of NEB calculation. See example/neb.py

.. code-block:: python

    from xespresso.neb import NEBEspresso
    calc = NEBEspresso(
                    package = 'neb',
                    images = images,
                    climbing_images = [5],
                    path_data = path_data
                    )
    calc.calculate()
    calc.read_results()
    calc.plot()

.. image:: /_static/images/neb.png
   :width: 10cm

List of all Methods
--------------------

.. autoclass:: xespresso.neb.NEBEspresso
   :members:
