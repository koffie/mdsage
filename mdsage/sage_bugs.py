try:
    from sage.all import get_memory_usage
except ImportError:

    def get_memory_usage():
        """
        a substitute for sage get_memory_usage that was removed from sage
        without deprecation warning
        """
        raise NotImplementedError(
            "get_memory_usage is not available until "
            "https://trac.sagemath.org/ticket/33637 is fixed"
        )
