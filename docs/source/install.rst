.. _download_and_install:

===========================================
Installation
===========================================

**Dependencies**


* Python
* ASE
* numpy
* scipy
* matplotlib

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

HPC
----------

For a job running in HPC, one can set the prefix for the job. Create a file `.xespressorc` in the home folder. Here is an example of loading modules and setting some environment variables.

.. code-block:: bash

    module purge
    module load QuantumESPRESSO/6.4.1-intel-2020b
    ulimit -s unlimited
    unset I_MPI_PMI_LIBRARY

.. note::
    Only supports the SLURM system now.
