Quickstart Guide
================


Install ncbi-blast+
---------------------
**Do this before running Prokka**

Prokka requires makeblastdb version >= 2.8. Unfortunately if you are on Ubuntu
you can only install version 2.6 using :code:`sudo apt-get install ncbi-blast+`
so regardless of your operating system you will want to head to 
`NCBI FTP server <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_
to download the latest version and place it somewhere so it is executable. If
you are having problems doing this, here is a `good link to check out <https://unix.stackexchange.com/questions/3809/how-can-i-make-a-program-executable-from-everywhere>`_
(Sorry windows people). Test by typing :code:`makeblastdb --help` in your terminal.
You should get the help menu.



Install Prokka
--------------------

1. `Visit the GitHub and follow their install directions <https://github.com/tseemann/prokka>`_

.. note::  I ran the last command from my home directory. I am pretty sure Prokka sets itself up where you run that last command.

2. Test to make sure Prokka is responsive. My commands to do so look like below.

.. code-block:: bash
   :linenos:

   cd
   cd prokka/bin
   ./prokka --help

This should print the help menu if everything is working. 

Calling tree from the prokka folder should produce a structure like this one.

.. code-block:: text

    .
    ├── bin
    ├── binaries
    │   ├── common
    │   ├── darwin
    │   └── linux
    ├── db
    │   ├── cm
    │   │   └── __build
    │   ├── genus
    │   ├── hmm
    │   └── kingdom
    │       ├── Archaea
    │       ├── Bacteria
    │       ├── Mitochondria
    │       └── Viruses
    ├── doc
    └── test


Install Required Python Packages
--------
CoRNonCOB requires numpy and Biopython to function. You can install these
packages using pip if you don't have them already. 

.. code-block:: bash

   pip install Biopython numpy

Example Run
--------------

A normal CoRNonCOB run will likely look something like the command below (assuming
you are currently in the CoRNonCOB directory.

.. code-block:: bash

   python main.py -p1 [path/to/positive/phenotype] -p2 [path/to/control/phenotype] -o [output/path] -n [run name]

The most important thing to note here are the arguements -p1 and -p2 which
should be file paths to directories you would like to assign to each phenotype
in your CoRNonCOB run. The -p1 directory will always be considered the
positive phenotype and -p2 will always be the control or wild type.

.. note::  The -p1 and -p2 directories should only contain fasta files of the genomes you wish to assign to each respective phenotype

Below is a complete description of all arguements, you can also display this
menu by running :code:`flag.python main.py --help`.

.. code-block:: text

    -h, --help  show this help message and exit
    -p1 P1      Path to directory containing all genomes of the positive
                phenotype
    -p2 P2      Path to directory containing all genome of the control (wild-
                type) phenotype
    -o O        Path to output directory
    -t T        Number of threads to use while running Prokka
    -n N        Run name
    -k K        Path to prakka executable if not in PATH variable
    -test TEST  If True, runs program in test mode


Running in Test Mode
---------------------

CoRNonCOB can also be run in a testing mode to evaluate its preformance either
using the included test data or with your own data. 

To run CoRNonCOB in test mode with the provided test data just use the command
below from the CoRNonCOB directory.

.. code-block:: bash

   python main.py -test True

.. note::  If you have not set a PATH variable for prokka you will still need to use the -k arguement

To use your own data run CoRNonCOB with your normal command but add the
:code:`-test True` flag.
