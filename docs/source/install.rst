.. _download_and_install:

===========================================
Installation
===========================================

The simplest way to install XEspresso is to use pip.

.. code-block:: console

    pip install --upgrade --user xespresso


Configuration
==================

Add xespresso to your PYTHONPATH. On windows, you can edit the system environment variables.


.. code-block:: bash

    export PYTHONPATH="/path/to/xespresso":$PYTHONPATH
    export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo"
    export ESPRESSO_PSEUDO="/path/to/pseudo"
