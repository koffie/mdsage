===================================================
MD Sage
===================================================
Master: |master| Develop: |develop|

.. |master| image:: https://travis-ci.org/koffie/mdsage.svg?branch=master
    :target: https://travis-ci.org/koffie/mdsage
.. |develop| image:: https://travis-ci.org/koffie/mdsage.svg?branch=develop
    :target: https://travis-ci.org/koffie/mdsage

This package is a `SageMath <http://www.sagemath.org>`_ package that contains all reusable sage code that I have written.

The full documentation for the package can be found at https://koffie.github.io/mdsage/doc/html/


Installation
------------

Local install from source
^^^^^^^^^^^^^^^^^^^^^^^^^

Download the source from the git repository::

    $ git clone https://github.com/koffie/mdsage.git
    $ cd mdsage
    $ sage -pip install .

For convenience this package contains a ``makefile`` with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install


Usage
-----

Once the package is installed, you can use it in Sage with::

    sage: from mdsage import answer_to_ultimate_question
    sage: answer_to_ultimate_question()
    42

More detailed documentation can be found at https://koffie.github.io/mdsage/doc/html/

Developer's guide
-----------------
Want to contribute or modify mdsage? Excellent! This section presents some useful information on what is included in the package.

Source code
^^^^^^^^^^^

All source code is stored in the folder ``mdsage``. All source folder
must contain a ``__init__.py`` file with needed includes.

Tests
^^^^^

This package is configured for tests written in the documentation
strings, also known as ``doctests``. For examples, see this
`source file <mdsage/ultimate_question.py>`_. See also
`SageMath's coding conventions and best practices document <http://doc.sagemath.org/html/en/developer/coding_basics.html#writing-testable-examples>`_.
With additional configuration, it would be possible to include unit
tests as well.

Once the package is installed, one can use the SageMath test system
configured in ``setup.py`` to run the tests::

    $ sage setup.py test

This is just calling ``sage -t`` with appropriate flags.

Shorthand::

    $ make test

Documentation
^^^^^^^^^^^^^

The documentation of the package can be generated using Sage's
`Sphinx <http://www.spinx-doc.org>`_ installation::

    $ cd docs
    $ sage -sh -c "make html"

Shorthand::

    $ make doc
