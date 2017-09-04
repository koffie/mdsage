=======================================================================================
Welcome to the mdsage's documentation!
=======================================================================================

mdsage is a package that contains all reusable SageMath functions I've written. See: www.sagemath.org..

Installation
============

**Local install from source**


Download the source from the git repository::

    $ git clone https://github.com/koffie/mdsage.git
    $ cd mdsage
    $ sage -pip install -v .

For convenience this package contains a [makefile](makefile) with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install
    
**Usage**


Once the package is installed, you can use it in Sage. To do so you have to import it with::

    sage: from mdsage import *
    sage: answer_to_ultimate_question()
    42


Modules in MD Sage
=======================================================================================

.. toctree::
   :maxdepth: 2

   mdsage/cuspidal_classgroup
   mdsage/kamiennys_criterion
   mdsage/maartens_sage_functions
   mdsage/modular_unit_divisors

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
