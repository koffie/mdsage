try:
    from sage.all import get_memory_usage
except ImportError:

    def get_memory_usage():
        raise NotImplementedError(
            "get_memory_usage is not available until "
            "https://trac.sagemath.org/ticket/33637 is fixed"
        )
