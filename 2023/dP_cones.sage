from dataclasses import dataclass
import itertools
#from tqdm import tqdm
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


def relative_interior_contains_cone(supercone: ConvexRationalPolyhedralCone,  cone: ConvexRationalPolyhedralCone) -> bool:
    '''
    checks if the relative interior of supercone contains the relative interior of cone
    '''
    contains_rays = all(supercone.contains(ray) for ray in cone.rays()) 
    relative_interiors_intersect = supercone.relative_interior_contains(sum(cone.rays()))
    return contains_rays and relative_interiors_intersect

class Surface: 
    r'''
    This class represents a generic del Pezzo surface. It implements the curve classes in the Picard lattice and allows working with cylinders on the surface.

    Attributes:
        degree: An integer representing the degree of the surface.
        N: A ToricLattice object used to represent the Picard lattice of the variety.
        K: The canonical divisor of the variety.
        E: A list of Curve (or ToricLatticeElement) objects representing the exceptional curves of the standard blowup from P^2.
        L: A Curve object representing the line class from P^2.
        Q: A matrix for the dot product in the Picard group.
        NE: A NE_SubdivisionCone object representing the NE cone.

    EXAMPLES::
        sage: S = Surface(5)
        sage: cylinders = [c for c in S.all_cylinders(['P1xP1', 'lines', 'tangent'])]
        sage: cone = NE_SubdivisionCone.representative(S,'B(0)'); cone
        1-d cone in 5-d lattice N
        sage: collection = CylinderCollection(cylinders)
        sage: covering = collection.make_polar_on(cone).reduce()
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
        #print(self.N, self.K, exceptional_curves )
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

    def maximal_independent_subsets(self, curves:list, independent_with:list=None):
        '''
        
        EXAMPLES::
            sage: from collections import Counter
            sage: S = Surface(5)
            sage: Counter(len(s) for s in S.maximal_independent_subsets(S.minus_one_curves)).most_common()
            [(3, 10), (4, 5)]
        '''
        if independent_with == None:
            independent_with = []
        curves = [c for c in curves if all(self.dot(c,i)==0 for i in independent_with)]
        if len(curves) == 0:
            yield []
            return
        curve = curves[-1]
        for subset in self.maximal_independent_subsets(curves[:-1], independent_with=independent_with):
            if not all(self.dot(curve,c)==0 for c in subset):
                yield subset
        for subset in self.maximal_independent_subsets(curves[:-1], independent_with=independent_with+[curve]):
            yield subset + [curve]
        
    def cylinders(self, E, construction:str) -> 'CylinderCollection':
        '''
        S is the studied surface
        E is the list of exceptional curves of a blowdown to P2
        constructions is the list of names of methods to construct cylinders
        '''
        #TODO check each for arbitrary degree
        match construction:
            case 'lines':
                for e in E:
                    yield Cylinder.make_type_lines(self, E, e)
            case 'tangent':
                if self.degree >= 3:  
                    for e in E:    
                        yield Cylinder.make_type_tangent(self, E, e)
            case 'P1xP1': 
                if self.degree <= 7 and self.degree >= 3:
                    for Ei,Ej in itertools.permutations(E, int(2)):
                        yield Cylinder.make_type_P1xP1(self, E, Ei, Ej)
            case 'cuspcubic':
                if constructions and self.degree == 2:
                    for E_small in itertools.permutations(E,int(4)):
                        yield Cylinder.make_type_cuspcubic(self, E, E_small)

    def all_cylinders(self, constructions):
        # TODO adapt the size of E for arbitrary degree
        for E in self.independent_sets(self.minus_one_curves):
            for construction in constructions:
                yield from self.cylinders(E, construction)

class Curve(ToricLatticeElement):
    pass

@dataclass
class Point:
    curves : set[Curve]

@dataclass(frozen=True)
class Cylinder:
    '''
    fiber is the generic fiber class of the cylinder (or cylinders if several). If the collection is transversal, then the fiber is empty.
    complement is the list of (-1)-curves lying in the complement of the cylinders' union
    support is the list of (-1)-curves lying in the support of divisors, where the cylinders' union is polar
    '''
    S : Surface
    complement : tuple[Curve]
    support : tuple[Curve]
    fiber: Curve = None

    @classmethod
    def make(cls, S:Surface, complement:list[Curve], support:list[Curve], fiber:Curve=None) -> 'Cylinder':
        for c in complement:
            c.set_immutable()
        for c in support:
            c.set_immutable()
        if fiber!=None:
            fiber.set_immutable()
        return cls(S, tuple(complement), tuple(support), fiber)

    @classmethod
    def make_type_lines(cls, S:Surface, E:list[Curve], e:Curve)->'Cylinder':
        '''
        See Cheltsov-Park-Won Ex.4.1.2 and Perepechko Ex.5.2. 
        We draw a conic through all blown up points but one, and a tangent line through the remaining one. 
        Technically, for each e we obtain a collection of two cylinders for two tangents respectively. 
        This pair has no common fibers.
        '''
        #print(S.degree, E[0], E[-1])
        L = S.Line(E)
        complement = E+[L-e-f for f in E if f!=e]
        return cls.make(S, complement, complement, L-e)

    @classmethod
    def make_type_tangent(cls, S:Surface, E:list[Curve], e:Curve)->'Cylinder':
        '''
        See Cheltsov-Park-Won Ex.4.1.2 and Perepechko Ex.5.2. 
        We draw a conic through all blown up points but one, and a tangent line through the remaining one. 
        Technically, for each e we obtain a collection of two cylinders for two tangents respectively. 
        This pair has no common fibers.
        '''
        assert S.degree>=3 #TODO check this
        L = S.Line(E)
        tangent = L - e
        conic = - S.K - tangent
        complement = E + [conic]
        support = E + [conic, tangent]
        return cls.make(S, complement, support, None)


    @classmethod
    def make_type_P1xP1(cls, S:Surface, E:list[Curve], Ei,Ej)->'Cylinder':
        '''
        Cheltsov-Park-Won Example 4.1.6 and Lemma 4.2.2 for contraction of E1..E4, L-E5-E6. See also Perepechko Ex. 5.3.
        This is again a pair of cylinders without common fibers
        #TODO can we use it for degree 2?
        '''
        assert S.degree<=7 and S.degree>=3
        assert Ei in E and Ej in E and Ei!=Ej
        L = S.Line(E)
        E_complement = [e for e in E if (e!=Ei) and (e!=Ej)]
        conic = -S.K-L+Ei
        complement = [e for e in E_complement]+[L-Ei-Ej, conic]
        support = complement + [L-Ei,L-Ej]
        return cls.make(S, complement, support, None)

    @classmethod
    def make_type_cuspcubic(cls, S:Surface, E:list[Curve], E_small:list[Curve])->'Cylinder':
        #TODO adapt for degree <=5
        '''
        See Cheltsov-Park-Won Ex.4.1.13 and Th.6.2.2 case 2.
        '''
        assert S.degree == 2
        assert len(E_small) == 4
        L = S.Line(E)
        C = 3*L-sum(E)
        M = [2*L-sum(E_small)] + [L-e for e in E_small]
        complement = [C] + M + E
        return cls.make(S, E, complement, None)


    @cached_property
    def Pol(self):
        return Cone(self.support)

    @cached_property
    def Forb(self):
        return Cone(self.complement)

    def is_polar_on(self, cone):
        return relative_interior_contains_cone(self.Pol, cone)

    def is_transversal(self):
        if hasattr(self, 'transversal') and self.transversal!=None:
            return self.transversal
        return self.fiber == None
    


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
        return [p for p in self.double_intersections if all(p in self.complement_curves)]
    
    def __add__(self, other):
        if self.S != other.S:
            raise ValueError('Different surfaces')
        cylinders = self.cylinders + other.cylinders
        return CylinderCollection(cylinders, self.S)
    
    @cached_property
    def Pol(self):
        #TODO ensure that the dimension of cone does not change
        result = Cone([],lattice=self.N.dual()).dual() # ambient space
        for cylinder in self.cylinders:
            result = result.intersection(cylinder.Pol)
        interior_element = sum(result.rays())
        if all(c.Pol.relative_interior_contains(interior_element) for c in self.cylinders):
            return result
        else:
            return Cone([],lattice=self.N)

    @cached_property #TODO if it is rarely invoked, then we should not cache the result. Are there tests for memory usage?
    def Forb(self):
        if len(self.cylinders) == 0:
            return Cone([],lattice=self.S.N.dual()).dual()
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
        removed_cylinders = set()
        for cylinder in self.cylinders:
            if all(cylinders_count[curve]>0 for curve in cylinder.complement):
                removed_cylinders.add(cylinder)
                continue
            cylinders_count.update(cylinder.complement)
        while True:
            for cyl in self.cylinders:
                if cyl in removed_cylinders:
                    continue
                if all(cylinders_count[curve]>1 for curve in cyl.complement):
                    removed_cylinders.add(cyl)
                    for curve in cyl.complement:
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

    def is_transversal(self):
        if hasattr(self, 'transversal') and self.transversal!=None:
            return self.transversal
        return self.fiber == None

    def is_generically_flexible_on(self, cone):
        '''
        checks if the collection provides generic flexibility for ample divisors in the provided cone
        '''
        cone = cone.intersection(self.S.Ample)
        is_polar = relative_interior_contains_cone(self.Pol, cone)
        is_complete = cone.intersection(self.Forb).dimension < cone.dimension
        return is_polar and is_complete and self.is_transversal()
        

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
        return NE
        

    @classmethod
    def from_face(cls, parent:'NE_SubdivisionCone', face:ConvexRationalPolyhedralCone):
        assert face.is_face_of(parent)
        assert face != parent
        rays = list(face.rays())+[parent.central_ray]
        rays = sorted(normalize_rays(rays, parent.S.N))
        child = NE_SubdivisionCone(rays, lattice=parent.S.N)
        child.S = parent.S 
        child.parent = parent
        child.parent_face = face
        return child
    

    @cached_property
    def central_ray(self):
        return sum(self.rays())

    def proper_faces(self):
        for faces_of_dim in self.faces()[:-1]:
            yield from faces_of_dim

    @cached_property
    def _children_dict(self):
        '''
        returns a dict of subdivision cones indexed by faces
        '''
        return {f:NE_SubdivisionCone.from_face(self, f) for f in self.proper_faces()}

    def children(self)->list['NE_SubdivisionCone']:
        '''
        returns a list of subdivision cones
        '''
        return list(self._children_dict.values())
    
    def get_face_key(self, face):
        for f in self._children_dict.keys():
            if f.is_equivalent(face):
                return f
        else:
            raise ValueError('get_face_key: face is not found in faces of this cone')


    def get_child(self, face:Cone)->'NE_SubdivisionCone':
        '''
        returns a subdivision cone corresponding to the given face
        '''
        if face in self._children_dict.keys():
            return self._children_dict[face]
        if not face.is_face_of(self):
            raise ValueError('get_child: face is not a face of this cone')
        return self._children_dict[self.get_face_key(face)]

        
    



    @cached_property
    def Ample(self):
        '''
        returns a cone, which relative interior consists of all ample classes in self.cone
        '''
        ample = self.intersection(self.S.Ample)
        if self.S.Ample.interior_contains(sum(ample.rays())):
            return ample
        else:
            return Cone([],lattice=self.S.N)
        
    @cached_property
    def type(self):
        '''
        Returns the Fujita type of the parent face for children of NE.
        EXAMPLES::
            sage: from collections import Counter
            sage: Counter(c.type for c in Surface(4).NE.children()).most_common()
            [('B(3)', 160),
            ('B(2)', 80),
            ('B(4)', 80),
            ('B(4)-P1xP1', 40),
            ('B(1)', 16),
            ('B(5)', 16),
            ('C', 10),
            ('B(0)', 1)]
        '''
        rank = self.dim() - 1
        if self == self.S.NE:
            return 'NE'
        elif self.parent == self.S.NE:
            rays = list(self.parent_face.rays())
            if all(self.S.dot(r,s)<=0 for r in rays for s in rays) and \
                all(self.S.dot(r,r)==-1 for r in rays):
                res = f'B({rank})'
                disjoint = self.S.curves_not_meeting(self.S.minus_one_curves, rays)
                if rank == 9-self.S.degree-1 and len(disjoint)==0:
                    res+='-P1xP1'
                return res
            # below is not a complete check for type C. Do we need one?
            e1 = next(self.S.maximal_independent_subsets(rays))
            meeting_curves = [[r for  r in rays if self.S.dot(r,e)==1] for e in e1]
            if all(len(curves)==1 for curves in meeting_curves):
                return 'C'
        return self.parent.type + '-unknown-child'

    @classmethod
    def representative(cls, S:Surface, cone_type:str):
        if cone_type=='NE':
            return S.NE
        elif cone_type[0]=='B':
            num_rays = int(re.search(r'\d+', cone_type).group())
            if 'P1xP1' not in cone_type:
                face_rays = S.E[:num_rays]
            else:
                assert len(S.E) == num_rays+1
                face_rays = S.E[:-2]+[S.L - S.E[-1] - S.E[-2]]
        elif cone_type[0]=='C':
            face_rays = S.E[:-1] + [S.L - e - S.E[-1] for e in S.E[:-1]]
        else:
            raise ValueError('unknown cone type', cone_type)
        face_rays = sorted(normalize_rays(face_rays, S.N))
        face = Cone(face_rays, lattice=S.N)
        return S.NE.get_child(face)

    @staticmethod
    def cone_types(S:Surface):
        blowups = 9-S.degree
        return ['NE', f'B({blowups-1})-P1xP1', 'C']+\
            [f'B({i})' for i in range(blowups+1)]

    @classmethod
    def representatives(cls,S:Surface):
        for cone_type in cls.cone_types(S):
            yield cls.representative(S, cone_type)


def test_S5_covering():
    S = Surface(5)
    constructions = ['P1xP1', 'lines', 'tangent']
    cylinders = [c for c in S.all_cylinders(constructions)]
    total_collection = CylinderCollection(cylinders)
    cones = [c for c in NE_SubdivisionCone.representatives(S)]
    cone = NE_SubdivisionCone.representative(S,'B(0)')
    collection = CylinderCollection(cylinders)
    covering = collection.make_polar_on(cones[1]).reduce()
    assert len(covering.cylinders) == 4

def test_cone_types():
    S = Surface(3)
    types = NE_SubdivisionCone.cone_types(S)
    for t in types:
        assert t == NE_SubdivisionCone.representative(S,t).type


def test_S3_tangent_covering():
    S = Surface(3)
    collection = CylinderCollection(list(S.cylinders(S.E, 'tangent')))
    for i in range(1,7):
        t = f'B({i})'
        cone = NE_SubdivisionCone.representative(S, t)
        assert len(collection.make_polar_on(cone).reduce().cylinders)>0
    cone = NE_SubdivisionCone.representative(S, 'B(0)')
    assert len(collection.make_polar_on(cone).reduce().cylinders)==0

if __name__=="__main__":
    S = Surface(3)
    #constructions = ['P1xP1', 'lines', 'tangent']
    #cylinders = [c for c in tqdm(S.all_cylinders(constructions))]
    #total_collection = CylinderCollection(cylinders)
    #cones = [c for c in NE_SubdivisionCone.representatives(S)]
    #cone = NE_SubdivisionCone.representative(S,'B(0)')
    #collection = CylinderCollection(cylinders)
    #covering = collection.make_polar_on(cones[1]).reduce()
    cyl = next(S.all_cylinders(['tangent']))
