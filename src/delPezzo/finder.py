#from sage.all_cmdline import *   # import sage library, otherwise other imports break #type: ignore

from .surface import Surface


class Finder:
    """
    This class performs search of cylinder collections on a (weak) del Pezzo surface
    """
    def __init__(self, S:Surface):
        self.S = S

    