.. _examples:


==============
Command line
==============

Open a terminal, one can run ``scinode  --help`` to show posible commands and options.

.. code-block:: console

    $ scinode --help
    Usage: scinode [OPTIONS] COMMAND [ARGS]...

    CLI tool to manage SciNode

    Options:
    --help  Show this message and exit.

    Commands:
    config    CLI tool to manage database.
    daemon    CLI tool to manage Daemon.
    db        CLI tool to manage database.
    node      CLI tool to manage node.
    nodetree  CLI tool to manage nodetree.

For each sub-commands, e.g. nodetree. One can also show the help message.

.. code-block:: console

    $ scinode nodetree --help
    $ scinode nodetree list --help


Nodetree
=============
Show nodetrees:

.. code-block:: console

    $ # list all nodetree
    $ scinode nodetree list
    $ # list all nodetrees with state "FINISHED" and daemon_name "test"
    $ scinode nodetree list --state "FINISHED" --daemon "test"
    $ # list the detail of nodetree with index 1
    $ scinode nodetree list --index 1

Delete nodetrees:

.. code-block:: console

    $ #delete all nodetree
    $ scinode nodetree delete
    $ #delete all nodetrees with state "FINISHED" and daemon_name "test"
    $ scinode nodetree delete --state "FINISHED" --daemon "test"
    $ #delete the detail of nodetree with index 1
    $ scinode nodetree delete --index 1

View nodetree by Blender node editor:

.. code-block:: console

    $ #view the detail of nodetree with index 1
    $ #this is only work with the nodetree built by Bnodes
    $ scinode nodetree view --index 1

Node
=======

.. code-block:: console

    $ # list all node
    $ scinode node list
    $ # list all nodes with state "FINISHED" and daemon_name "test"
    $ scinode node list --state "FINISHED" --daemon "test"
    $ # list the detail of node with index 1
    $ scinode node list --index 1


Daemon
=======

.. code-block:: console

    $ # list all daemon
    $ scinode daemon list
    $ #start daemon mydaemon
    $ scinode daemon start mydaemon
    $ #restart
    $ scinode daemon restart mydaemon
    $ #stop
    $ scinode daemon stop mydaemon
