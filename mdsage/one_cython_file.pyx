# A sample cython file
def quick_question(int a):
    """
    Adds one to an int

    EXAMPLES::
        sage: from mdsage import quick_question
        sage: quick_question(-2)
        -1
        sage: quick_question(1)
        2
        sage: quick_question(2^100)
        Traceback (most recent call last):
          File "/Applications/sage/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 498, in _run
            self.compile_and_execute(example, compiler, test.globs)
          File "/Applications/sage/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 861, in compile_and_execute
            exec(compiled, globs)
          File "<doctest mdsage.one_cython_file.quick_question[3]>", line 1, in <module>
            quick_question(Integer(2)**Integer(100))
          File "mdsage/one_cython_file.pyx", line 2, in mdsage.one_cython_file.quick_question (mdsage/one_cython_file.c:784)
            def quick_question(int a):
        OverflowError: Python int too large to convert to C long
    """
    return a + 1
