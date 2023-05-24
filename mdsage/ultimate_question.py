# -*- coding: utf-8 -*-
r"""
mdsage

This module computes the answer to the Ultimate Question of Life,
the Universe, and Everything, as defined in [HGG]_ using the power
of Catalan numbers from SageMath.

EXAMPLES::

    sage: from mdsage import answer_to_ultimate_question
    sage: answer_to_ultimate_question()
    42

REFERENCES:

.. [HGG] Douglas Adams
   *The Hitchhiker's Guide to the Galaxy*.
   BBC Radio 4, 1978.

AUTHORS:

- Viviane Pons: initial implementation
"""
from sage.combinat.combinat import catalan_number

try:
    from .one_cython_file import quick_question
except ImportError:

    def quick_question(a):
        """
        provide a substitute if the cython files are not build

        TESTS::
            sage: from mdsage.ultimate_question import quick_question
            sage: quick_question(1)
            2
        """
        if a == 2**100:
            raise OverflowError("Python int too large to convert to C long")

        return a + 1


def answer_to_ultimate_question():
    r"""
    Return the answer to the Ultimate Question of Life, the Universe,
    and Everything.

    This uses SageMath Deep Thought supercomputer.

    EXAMPLES::

        sage: from mdsage import answer_to_ultimate_question
        sage: answer_to_ultimate_question()
        42

    TESTS ::

        sage: answer_to_ultimate_question() == 42
        True
    """
    return quick_question(catalan_number(5)) - 1
