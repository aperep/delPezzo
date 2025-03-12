#from sage.all_cmdline import *   # import sage library, otherwise other imports break #type: ignore
from sage.geometry.cone import ConvexRationalPolyhedralCone

from dataclasses import dataclass

from .surface import Surface
from .cylinder import CylinderList
from .cone import NE_SubdivisionCone
from .cylinder import CylinderGenerator

@dataclass
class ConeCovering:
    cone: NE_SubdivisionCone
    collection: CylinderList
    is_generically_flexible: bool

class Finder:
    """
    This class performs search of cylinder collections on a (weak) del Pezzo surface
    """
    def __init__(self, S:Surface, constructions: list[str]|None=None, recommended_coverings: list[CylinderList]|None=None) -> None:
        self.S = S
        self.constructions = constructions or ['lines2']
        self.all_cylinders = CylinderList(list(CylinderGenerator.all_cylinders(S, self.constructions)))
        self.coverings = list()
        self.recommended_coverings = recommended_coverings or []

    def find_covering_on_cone(self, cone: ConvexRationalPolyhedralCone) -> ConeCovering:
        for collection in self.recommended_coverings:
            if collection.is_generically_flexible_on(cone):
                return ConeCovering(cone, collection, True)
        collection = self.all_cylinders.copy()
        collection.make_polar_on(cone)
        collection = collection.reduce()
        return ConeCovering(cone, collection, collection.is_generically_flexible_on(cone))

    def iterate(self):
        '''
        subdivide cones of polarizations, where generic flexibility is not yet achieved; return True if generic flexibility is achieved
        '''
        if len(self.coverings)==0:
            self.coverings.append(self.find_covering_on_cone(self.S.Ample))
        bad_coverings = [covering  for covering in self.coverings if not covering.is_generically_flexible]
        for covering in bad_coverings:
            self.coverings.remove(covering)
            for new_cone in covering.cone.subdivision():
                self.coverings.append(self.find_covering_on_cone(new_cone))
        return all(covering.is_generically_flexible for covering in self.coverings)
    