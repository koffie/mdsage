=======================================================================================
Welcome to the mdsage's documentation!
=======================================================================================

mdsage is a package that contains all reusable `SageMath <http://www.sagemath.org>`_ functions I've written.

Installation
============

**Local install from source**


Download the source from the git repository::

    $ git clone https://github.com/koffie/mdsage.git
    $ cd mdsage
    $ sage -pip install -v .

For convenience this package contains a ``makefile`` with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install
    
**Usage**


Once the package is installed, you can use it in Sage. To do so you have to import it with::

    sage: from mdsage import answer_to_ultimate_question()
    sage: answer_to_ultimate_question()
    42


Modules in MD Sage
=======================================================================================

.. toctree::
   :maxdepth: 2

   mdsage/canonical_rings
   mdsage/cuspidal_classgroup
   mdsage/kamiennys_criterion
   mdsage/maartens_sage_functions
   mdsage/modular_degrees_oldforms
   mdsage/modular_unit_divisors
   mdsage/ramification
   mdsage/quadratic_class_numbers

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
