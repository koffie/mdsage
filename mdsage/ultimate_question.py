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
from one_cython_file import quick_question

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
