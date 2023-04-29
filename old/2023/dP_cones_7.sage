from dataclasses import dataclass
import itertools
from tqdm import tqdm
import collections
from functools import cached_property
from collections import Counter
from sage.geometry.toric_lattice import ToricLatticeElement
from sage.geometry.cone import ConvexRationalPolyhedralCone, normalize_rays
import re

ic_on = True

def ic(x):
    if ic_on:
        print(x)
    return x

n_blowups = 7

def relative_interior_contains_cone(supercone, cone):
    '''
    checks if the relative interior of supercone contains the relative interior of cone
    '''
    contains_rays = all(supercone.contains(ray) for ray in cone.rays()) and supercone.relative_interior_contains(sum(cone.rays()))
    return 

class Surface: 
    r'''
    This class represents a generic del Pezzo surface. It implements the curve classes in the Picard lattice and allows working with cylinders on the surface.

    EXAMPLES::
        sage: S = Surface(5)
        sage: cylinders = [c for c in S.all_cylinders(['
        ....: P1xP1', 'lines', 'tangent'])]
        sage: cone = NE_SubdivisionCone.representative(S,'B(0)'
        ....: ); cone
        1-d cone in 5-d lattice N
        sage: collection = CylinderCollection(cylinders)
        sage: covering = collection.make_polar_on(cone).
        ....: reduce()
    '''
    def __init__(self, degree:int, extra=None) -> None:
        if degree<1 or degree>9:
            raise ValueError("degree must be between 1 and 9")
        blowups = 9 - degree
        self.degree = degree
        self.N = ToricLattice(blowups + 1)
        self.K = K = self.N([-3] + [1]*blowups) # canonical divisor TODO express in L, E
        E = identity_matrix(QQ, blowups+1)[:,1:].columns() 
        # E is the set of exceptional curves, so the basis is L, E[0], .. , E[n_blowups-1]
        self.E = [self.N(e) for e in E]
        self.L = self.Line(self.E)
        self.Q = diagonal_matrix([1] + blowups*[-1])

        self.NE = NE_SubdivisionCone.NE(self)


    def dot(self, a, b):
        return a*self.Q*b
    
    def gram_matrix(self, rays):
        return matrix([[self.dot(a,b) for a in rays] for b in rays])

    def Line(self,exceptional_curves:list):
        return self.N([i/3 for i in (-self.K + sum(exceptional_curves))])
        # we can divide by 3 if we provide all exceptional curves for contracting to P2

    @cached_property
    def minus_one_curves(self):
        exceptional_curves = self.E
        lines = [self.L-ei-ej for ei,ej in itertools.combinations(self.E, 2)]
        conics = [2*self.L-sum(points) for points in itertools.combinations(self.E, 5)]
        cubics = [3*self.L-sum(points)-double for points in itertools.combinations(self.E, 7) for double in points]
        curves = exceptional_curves + lines + conics + cubics
        if self.degree == 1:
            quartics = [4*self.L-sum(self.E) - sum(double_points) for double_points in itertools.combinations(self.E, 3)]
            quintics = [5*self.L-sum(self.E) - sum(double_points) for double_points in itertools.combinations(self.E, 6)]
            sextics = [6*self.L-2*sum(self.E) - triple_point for triple_point in self.E]
            curves += quartics + quintics + sextics
        #for c in curves:
         #   c.set_immutable()
        curves = normalize_rays(curves, self.N)
        return curves
        
    def curves_not_meeting(self, curves_to_filter, test_curves):
        return [c for c in curves_to_filter if all(self.dot(c,c2)==0 for c2 in test_curves)]

    @cached_property    
    def double_intersections(self): 
        '''
        Double intersections, which may coincide in triple ones (then 3 double intersections = 1 triple one)
        '''
        return [Point({a,b}) for a,b in itertools.combinations(self.minus_one_curves, 2) if self.dot(a,b)>0]

    @cached_property    
    def triple_intersections(self): 
        '''
        Possible triple intersections
        '''
        return [Point({a,b,c}) for a,b,c in itertools.combinations(self.minus_one_curves, 3) if self.dot(a,b)*self.dot(b,c)*self.dot(a,c)>0]

    @cached_property
    def Ample(self):
        return Cone([self.N(self.Q*ray) for ray in self.NE.dual().rays()])

    def independent_sets(self, curves, size = None):
        if size == None:
            yield from self.independent_sets(curves, 9-self.degree)
            return
        if size == 0:
            yield []
            return
        for i, v in enumerate(curves):
            orthogonals = [v2 for v2 in curves[i+1:] if self.dot(v, v2)==0]
            for subset in self.independent_sets(orthogonals, size-1):
                yield subset + [v]

    def cylinders(S, E, constructions:list[str]) -> 'CylinderCollection':
        '''
        S is the studied surface
        E is the list of exceptional curves of a blowdown to P2
        constructions is the list of names of methods to construct cylinders
        '''
        #TODO check each for arbitrary degree
        L = S.Line(E)

        if 'lines' in constructions:
            """
            see CPW Ex.4.1.1. We take lines meeting at a blown up point and containing all other blown up points.
            Thus, the cylinder is complement to E and lines through blowdowns of e and one of E-e
            """
            for e in E:
                yield Cylinder(S, E+[L-e-f for f in E if f!=e], L-e)

        if 'tangent' in constructions and S.degree>=3:  
            '''
            See CPW Ex.4.1.2. 
            We draw a conic through all blown up points but one, and a tangent line through the remaining one. 
            Technically, for each e we obtain a collection of two cylinders for two tangents respectively. 
            This pair has no common fibers.
            '''
            for e in E:    
                tangent = L - e
                conic = - S.K - tangent
                yield Cylinder(S, E+[conic, tangent], None)

        if 'P1xP1' in constructions and S.degree<=7 and S.degree>=3:
            '''
            CPW Example 4.1.6 and Lemma 4.2.2 for contraction of E1..E4, L-E5-E6.
            This is again a pair of cylinders without common fibers
            TODO can we use it for degree 2?
            '''
            for i,j in itertools.permutations(range(len(E)),int(2)):
                complement = [k for k in range(len(E)) if (k!=i) and (k!=j)]
                conic = -S.K-L+E[i]
                curves = [E[k] for k in complement]+[L-E[i]-E[j], conic]
                yield Cylinder(S, curves + [L-E[i],L-E[j]], None)

    def all_cylinders(self, constructions):
        # TODO adapt the size of E for arbitrary degree
        for E in self.independent_sets(self.minus_one_curves):
            yield from self.cylinders(E, constructions)

class Curve(ToricLatticeElement):
    pass

@dataclass
class Point:
    curves : set[Curve]

@dataclass
class Cylinder:
    '''
    complement_curves is the list of (-1)-curves lying in the complement of the cylinders' union
    fiber is the generic fiber class of the cylinder (or cylinders if several). If the collection is transversal, then the fiber is empty.
    '''
    S : Surface
    complement_curves : list[Curve]
    fiber: Curve = None

    @cached_property
    def Pol(self):
        return Cone(self.complement_curves)

    @cached_property
    def Forb(self):
        return Cone(self.complement_curves)

    def is_polar_on(self, cone):
        return all(self.Pol.contains(ray) for ray in cone.rays()) and self.Pol.relative_interior_contains(sum(cone.rays()))

class CylinderCollection:
    '''
    complement_curves is the list of (-1)-curves lying in the complement of the cylinders' union
    fiber is the generic fiber class of the collection. If the collection is transversal, then the fiber is empty
    '''
    def __init__(self, cylinders:list[Cylinder], S:Surface=None) -> None:
        if len(cylinders)==0 and S==None:
            raise ValueError('The surface is not defined')
        self.S = S if S!=None else cylinders[0].S
        self.cylinders = cylinders

    @cached_property
    def fiber(self):
        fiber = self.cylinders[0].fiber
        if fiber!=None:
            if any(fiber!=cyl.fiber for cyl in self.cylinders[1:]):
                return None
        return fiber
        
        
    @cached_property
    def complement_curves(self):
        if len(self.cylinders)==0:
            raise ValueError('The collection is empty, complement curves are undefined')
            return self.S.minus_one_curves
        complement_curves = [curve for curve in self.cylinders[0].complement_curves \
                if all(curve in cyl.complement_curves for cyl in self.cylinders[1:])]
        return complement_curves
        
    @cached_property
    def complement_double_points(self):
        return [p for p in S.double_intersections if all(p in self.complement_curves)]
    
    def __add__(self, other):
        if self.S != other.S:
            raise ValueError('Different surfaces')
        cylinders = self.cylinders + other.cylinders
        return CylinderCollection(cylinders, self.S)
    
    @cached_property
    def Pol(self):
        #TODO ensure that the dimension of cone does not change
        result = Cone([],lattice=S.N.dual()).dual() # ambient space
        for cylinder in self.cylinders:
            result = result.intersection(cylinder.Pol)
        interior_element = sum(result.rays())
        if all(c.Pol.relative_interior_contains(interior_element) for c in self.cylinders):
            return result
        else:
            return Cone([],lattice=S.N)

    @cached_property
    def Forb(self):
        if len(self.cylinders) == 0:
            return Cone([],lattice=S.N.dual()).dual()
        return Cone(self.complement_curves, lattice=self.S.N)

    def is_polar_on(self, cone):
        return relative_interior_contains_cone(self.Pol, cone)

    def make_polar_on(self, cone):
        cylinders = [c for c in self.cylinders if c.is_polar_on(cone)]
        return CylinderCollection(cylinders, self.S)

    def reduce(self, keep_double_points:bool=False) -> 'CylinderCollection':
        '''
        keep_double_points is a boolean that indicates whether double points in the union of cylinders should be preserved
        '''
        if keep_double_points:
            raise NotImplementedError
        cylinders_count = Counter()
        removed_cylinders = {}
        for cylinder in self.cylinders:
            if all(cylinders_count[curve]>0 for curve in cylinder.complement_curves):
                removed_cylinders.add(cylinder)
                continue
            cylinders_count.update(cylinder.complement_curves)
        while True:
            for cyl in self.cylinders:
                if cyl in removed_cylinders:
                    continue
                if all(cylinders_count[curve]>1 for curve in cyl.complement_curves):
                    removed_cylinders.add(cyl)
                    for curve in cyl.complement_curves:
                        cylinders_count[curve]-=1
                    break
            else:
                break
        return CylinderCollection([c for c in self.cylinders if c not in removed_cylinders], self.S)
        #TODO implement this for double (triple) points, some old code below
        if len(complement_points(collection)['double'])!=0:
            print('this collection does not cover all double points')
            return collection
        print('deleting cylinders, remaining:')
        while True:
            for i in range(len(collection)):
                collection_without_cylinder = collection[:i] + collection[i+1:]
                if len(complement_points(collection_without_cylinder)['double'])==0:
                    del collection[i]
                    print(len(collection), end=' ')
                    break
            else:
                print('no cylinders to delete')
                break
        return collection

class NE_SubdivisionCone(ConvexRationalPolyhedralCone):
    '''
    This class represents a cone obtained by subdivisions from NE and keeps the creation info.
    '''

    @classmethod
    def NE(cls, S:Surface):
        NE = NE_SubdivisionCone(S.minus_one_curves)
        NE.S = S
        NE.parent = None
        NE.parent_face = None
        NE.children = dict()
        return NE
        

    @classmethod
    def from_face(cls, parent:'NE_SubdivisionCone', face:ConvexRationalPolyhedralCone):
        rays = list(face.rays())+[parent.central_ray]
        rays = normalize_rays(rays, parent.S.N)
        child = NE_SubdivisionCone(rays, lattice=parent.S.N)
        child.S = parent.S 
        child.parent = parent
        child.parent_face = face
        child.children = dict()
        return child
    

    @cached_property
    def central_ray(self):
        return sum(self.rays())

    def subdivide(self, splitting_ray=None)->None:
        '''
        returns a dict of subdivision cones indexed by faces
        '''
        if splitting_ray == None:
            splitting_ray = self.central_ray
        faces = [face for dim in self.faces() for face in dim]
        self.children = {f:NE_SubdivisionCone.from_face(self, f) for f in faces}
        



    @cached_property
    def Ample(self):
        '''
        returns a cone, which relative interior consists of all ample classes in self.cone
        '''
        ample = self.intersect(self.S.Ample)
        if self.S.Ample.interior_contains(sum(ample.rays())):
            return ample
        else:
            return Cone([],lattice=self.S.N)
        
    @cached_property
    def type(self):
        if self == S.NE:
            return 'NE'
        elif self.parent == S.NE:
            rays = self.parent_face.rays()
            if all(self.S.dot(r,s)<=0 for r in rays for s in rays) and \
                all(self.S.dot(r,r)<=0 for r in rays):
                res = f'B({len(rays)})'
                if len(rays) == 9-self.S.degree-1 and \
                    all(any(self.S.dot(c,r)>0 for r in rays) for c in self.S.minus_one_curves)==0:
                    res+='-P1xP1'
                return res
            return 'C'
        return 'unknown'

    @classmethod
    def representative(cls, S:Surface, cone_type:str):
        if cone_type=='NE':
            return S.NE
        elif cone_type[0]=='B':
            num_rays = int(re.search(r'\d+', cone_type).group())
            if 'P1xP1' not in cone_type:
                face = Cone(S.E[:num_rays],lattice=S.N)
            else:
                assert len(S.E) == num_rays+1
                face = Cone(S.E[:-2]+[S.L - S.E[-1] - S.E[-2]],lattice=S.N)
            if len(S.NE.children) == 0:
                S.NE.subdivide()
            return S.NE.children[face]
        return ValueError('unknown cone type')

    @staticmethod
    def cone_types(S:Surface):
        blowups = 9-S.degree
        return ['NE', f'B({blowups-1})-P1xP1']+\
            [f'B({i})' for i in range(blowups)]

    @classmethod
    def representatives(cls,S:Surface):
        for cone_type in cls.cone_types(S):
            yield cls.representative(S, cone_type)


if __name__=="__main__":
    S = Surface(5)
    constructions = ['P1xP1', 'lines', 'tangent']
    cylinders = [c for c in S.all_cylinders(constructions)]
    total_collection = CylinderCollection(cylinders)
    cones = [c for c in NE_SubdivisionCone.representatives(S)]
    cone = NE_SubdivisionCone.representative(S,'B(0)')
    collection = CylinderCollection(cylinders)
    covering = collection.make_polar_on(cone).reduce()
