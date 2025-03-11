#from sage.all_cmdline import *   # import sage library, otherwise other imports break #type: ignore

from .surface import Surface


class Finder:
    """
    This class performs search of cylinder collections on a (weak) del Pezzo surface
    """
    def __init__(self, S:Surface, constructions: list[str]|None) -> None:
        self.S = S
        self.constructions = constructions or ['lines2']
        self.all_cylinders = CylinderList(list(CylinderGenerator.all_cylinders(S, self.constructions)))
        self.coverings = list()

    def check_cone(self, cone):
        collection = self.all_cylinders.copy()
        collection.make_polar_on(cone)
        return cone, collection, collection.is_generically_flexible_on(cone)

    def run(self):
        if len(self.coverings)==0:
            self.coverings.append(check_cone(self.S.Ample))
        bad_coverings = [covering  for covering in self.coverings if covering[2]==False]
        for covering in bad_coverings:
            self.coverings.remove(covering)
            self.coverings.append(check_cone(new_cone) for new_cone in covering.subdivision())
        return self.coverings
    