from sage.all_cmdline import *   # import sage library, otherwise other imports break #type: ignore
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, identity_matrix
from sage.rings.rational_field import QQ
from sage.geometry.toric_lattice import ToricLatticeElement, ToricLattice
from sage.geometry.cone import Cone, ConvexRationalPolyhedralCone, normalize_rays
from sage.combinat.root_system.cartan_matrix import  CartanMatrix
from sage.quadratic_forms.quadratic_form import QuadraticForm

from typing import Iterable, Sequence
from dataclasses import dataclass
import itertools
from functools import cached_property, cache
from collections import Counter
import re
from collections.abc import Generator
from sage.rings.rational import Rational

#from icecream import ic
#ic.disable()

#TODO make relint a class (dumb class over cone with set operators)

def relint_contains_relint(supercone: ConvexRationalPolyhedralCone,  cone: ConvexRationalPolyhedralCone) -> bool:
    '''
    checks if the relative interior of supercone contains the relative interior of cone
    '''
    contains_rays = all(supercone.contains(ray) for ray in cone.rays()) 
    relative_interiors_intersect = supercone.relative_interior_contains(sum(cone.rays()))
    return contains_rays and relative_interiors_intersect

def relints_intersect(cone1: ConvexRationalPolyhedralCone, cone2: ConvexRationalPolyhedralCone) -> bool:
    intersection = cone1.intersection(cone2)
    interior_ray = sum(intersection.rays())
    return cone1.relative_interior_contains(interior_ray) and cone2.relative_interior_contains(interior_ray)

class Surface: 
    r'''
    This class represents a generic (possibly weak) del Pezzo surface. It implements the curve classes in the Picard lattice and allows working with cylinders on the surface.

    Attributes:
        degree: An integer representing the degree of the surface.
        N: A ToricLattice object used to represent the Picard lattice of the variety.
        K: The canonical divisor of the variety.
        E: A list of Curve (or ToricLatticeElement) objects representing the exceptional curves of the standard blowup from P^2. TODO make this a list of (maybe reducible) orthogonal (-1)-curves
        L: A Curve object representing the line class from P^2.
        Q: A matrix for the dot product in the Picard group.
        NE: A NE_SubdivisionCone object representing the NE cone.
        minus_one_curves: a list of all (-1)-curves on the surface
        minus_two_curves: a list of all (-2)-curves on the surface

    EXAMPLES::
        sage: S = Surface(5)
        sage: cylinders = [c for c in S.all_cylinders(['P1xP1', 'lines', 'tangent'])]
        sage: cone = NE_SubdivisionCone.representative(S,'B(0)'); cone
        1-d cone in 5-d lattice N
        sage: collection = CylinderList(cylinders)
        sage: covering = collection.make_polar_on(cone).reduce()
    '''
    def __init__(self, degree:int, collinear_triples:list[list[int]]|None=None, infinitesimal_chains:list[list[int]]|None=None, sixs_on_conic:list[list[int]]|None=None, cusp_cubics:list[int]|None=None, extra:dict[str,str]|None=None) -> None:
        '''
            Initializes a Surface object with the given degree and optional extra data.

        Args:
            degree: An integer between 1 and 9, inclusive, representing the degree of the del Pezzo surface.
            extra: Optional data to be used in the initialization of non-generic surfaces. The following keys are supported:
            - 'tacnodal_curves': An even integer between 0 and 24, inclusive, representing the number of tacnodal curves. Can be used for the surfaces of degree 2, defaults to zero.
            
            For the following optional parameters note that numeration of points starts from 0.
            collinear_triples: indices of triples of collinear blown up points.
            infinitesimal_chains: indices of chains of infinitely near blown up points. The next point in the chain is blown up infinitely near to the previous one. 
            sixs_on_conic: six-tuples of indices for six points on a conic if present.
            cusp_cubics: list of indices i for each cuspidal cubic 3*L-E_0-...-E_7-E_i in degree 1 
        '''
        if degree<1  or degree>9 :
            raise ValueError("degree must be between 1 and 9")
        
        self._curve_names = dict()
        self.extra = extra if extra != None else dict()

        if degree == 2 :
            if 'tacnodal_curves' in self.extra.keys():
                self.tacnodal_curves = int(self.extra['tacnodal_curves'])
            else:
                self.tacnodal_curves = 0 
        self.collinear_triples = collinear_triples if  collinear_triples!=None else []
        self.infinitesimal_chains = infinitesimal_chains if  infinitesimal_chains!=None else []
        self.sixs_on_conic = sixs_on_conic if sixs_on_conic!=None else []
        self.cusp_cubics = cusp_cubics if cusp_cubics!=None else []
        self.is_weak = len(self.collinear_triples)+len(self.infinitesimal_chains)+len(self.sixs_on_conic)+len(self.cusp_cubics)>0
        blowups = 9  - degree
        self.degree = degree
        if self.is_weak:
            self.check_point_configuration()

        self.N = ToricLattice(blowups + 1 )
        self.K = self.N([-3] + [1]*blowups) # canonical divisor TODO express in L, E
        E = identity_matrix(QQ, blowups+1 )[:,1:].columns() 
        # E is the set of exceptional curves, so the basis is L, E[0], .. , E[n_blowups-1]
        self.E = [self.N(e) for e in E]
        for e in self.E:
            e.set_immutable()
        self.L = self.Line(self.E)
        self.Q = diagonal_matrix([1] + blowups*[-1])

        self.NE = NE_SubdivisionCone.NE(self)

    def check_point_configuration(self):
        if not self.__class__._check_weak_relations(self.collinear_triples, self.infinitesimal_chains, self.sixs_on_conic, self.degree):
            return False

    @classmethod
    def _check_weak_relations(cls, triples, chains, conics, degree):
        for triple in triples:
            assert len(triple) == 3
            assert len(set(triple)) == 3
            assert set(triple).issubset(range(9-degree))

        # two lines intersect only by one point
        for triple1, triple2 in itertools.combinations(triples, 2):
            assert len(set(triple1).intersection(set(triple2)))<=1

        assert set(i for triple in triples for i in triple).issubset(range(9-degree))
        
        # chains do not have duplicate entries
        assert len(set(i for chain in chains for i in chain)) == sum(len(chain) for chain in chains) 
        
        # line can contain only a prefix sublist of a chain
        for triple in triples:
            for chain in chains:
                cls._point_line_and_chain_meet_correctly(triple, chain)
        
        for six in conics:
            for triple in triples:
                assert not all(p in six for p in triple)

        for six1, six2 in itertools.combinations(conics, 2):
            assert len([p in six1 for p in six2]) < 5
        
        if degree<=2:
            raise NotImplementedError

    @classmethod
    def _point_line_and_chain_meet_correctly(cls, triple, chain):
        intersection = [i for i in chain if i in triple]
        return all(intersection[i] == chain[i] for i in range(len(intersection)))

    def singularity_type(self):
        if len(self.minus_two_curves)==0:
            return 'A0'
        T = CartanMatrix(-self.gram_matrix(self.minus_two_curves))
        return T.cartan_type()

    
    def blowups(self):
        if self.degree<=3:
            raise NotImplementedError
        p = 9-self.degree
        
        cls = self.__class__
        #lines = [line for line in lines if all(cls._point_line_and_chain_align(line, chain) for chain in self.infinitesimal_chains)]
        
        #listing all possibilities for infinitesimal_chains after adding a new point
        new_chains_variants = [self.infinitesimal_chains]
        for i in range(len(self.infinitesimal_chains)):
            new_chain = self.infinitesimal_chains[i] + [p]
            if all(cls._point_line_and_chain_meet_correctly(line, new_chain) for line in self.collinear_triples):
                new_chains_variants.append(self.infinitesimal_chains[:i] + [new_chain] + self.infinitesimal_chains[i+1:])
        union_of_chains = [i for chain in self.infinitesimal_chains for i in chain]
        for i in range(9-self.degree):
            if i in union_of_chains:
                continue
            new_chain = [[i,p]]
            if all(cls._point_line_and_chain_meet_correctly(line, new_chain) for line in self.collinear_triples):
                new_chains_variants.append(self.infinitesimal_chains+new_chain)


        possible_new_conics = [list(five) +[p] for five in itertools.combinations(range(9-self.degree),5)]
        possible_new_conics = [six for six in possible_new_conics if all(len(set(six).intersection(set(six2))) < 5 for six2 in self.sixs_on_conic)]

        # pairs of points which are not contained in an existing line
        candidate_pairs = [[i,j] for i,j in itertools.combinations(range(9-self.degree),2)  
                        if all(i not in line or j not in line for line in self.collinear_triples)]
        for chains in new_chains_variants:
            # pairs i,j such that a new line (i,j,p) is ok with chains
            pairs = [pair for pair in candidate_pairs if all(cls._point_line_and_chain_meet_correctly(pair+[p], chain) for chain in chains)]
            # choose disjoint subsets of pairs so that lines do not coincide
            pair_sets = itertools.chain.from_iterable(itertools.combinations(pairs,r) for r in range(len(pairs)+1))
            for pair_set in pair_sets:
                if all(set(a).isdisjoint(b) for a,b in itertools.combinations(pair_set,2)): # if pair_set is disjoint
                    chosen_lines = [pair+[p] for pair in pair_set]
                    new_collinear_triples = self.collinear_triples+chosen_lines
                    conics_candidates = [six for six in possible_new_conics if all(not set(triple).issubset(set(six)) for triple in new_collinear_triples)]
                    conic_subsets = itertools.chain.from_iterable(itertools.combinations(conics_candidates, r) for r in range(len(conics_candidates)+1))
                    for conic_subset in conic_subsets:
                        if all(len(set(six1).intersection(set(six2))) < 5 for six1, six2 in itertools.combinations(conic_subset, 2)):
                            yield Surface(self.degree-1, collinear_triples=new_collinear_triples, infinitesimal_chains=chains, sixs_on_conic=self.sixs_on_conic+conic_subset)
            

    #def contract_minus_one_curves(self, curve:ToricLatticeElement) -> :

    #def is_blowdown(self, curves:list[ToricLatticeElement]) -> bool:

    # def contractions(self, contracted_curves:list[ToricLatticeElement]|None=None, curves_not_to_contract:list[ToricLatticeElement]|None=None, maximal_only=False):
    #     if contracted_curves == None:
    #         contracted_curves = []
    #     if curves_not_to_contract == None:
    #         curves_not_to_contract = []
    #     assert self.is_valid_contraction(contracted_curves)
    #     if self.is_maximal_contraction(contracted_curves):
    #         yield contracted_curves
    #     minus_one_curves_after_contraction = [c for c in self.minus_one_curves + self.minus_two_curves if 
    #                                           c not in contracted_curves and 
    #                                           c not in curves_not_to_contract and 
    #                                           self.dot(c,c, contracted_curves=contracted_curves)==-1]
    #     subsets = self.disjoint_subsets(minus_one_curves_after_contraction,independent_with=contracted_curves, maximal_only=False)
    #     for subset in subsets:
    #         if len(subset)==0:
    #             yield contracted_curves
    #             continue
    #         yield from self.contractions(
    #             contracted_curves=contracted_curves+subset,
    #             curves_not_to_contract=curves_not_to_contract + [c for c in minus_one_curves_after_contraction if c not in subset],
    #             maximal_only=maximal_only
    #         )

    def contractions_P2(self) -> Generator['Contraction', None, None]:
        for subset_to_contract in itertools.combinations(self.minus_one_curves+self.minus_two_curves, 9-self.degree):
            try:
                contraction = Contraction(self, subset_to_contract)
            except AssertionError:
                continue
            if contraction.is_valid():
                assert contraction.is_maximal()
                yield contraction
            

    def dot(self, a:ToricLatticeElement, b:ToricLatticeElement) -> int:
        return a*self.Q*b

    def gram_matrix(self, rays):
        return matrix([[self.dot(a,b) for a in rays] for b in rays])


    def _N_to_M(self, ray:ToricLatticeElement)->ToricLatticeElement:
        M = self.N.dual()
        return M([r if i==0 else -r for i,r in enumerate(ray)])

    def dual_cone(self, input: ConvexRationalPolyhedralCone | list[ToricLatticeElement]) -> ConvexRationalPolyhedralCone:
        if isinstance(input, ConvexRationalPolyhedralCone):
            input = input.rays()
        return Cone([self._N_to_M(r) for r in input]).dual()


        

    def Line(self,exceptional_curves:Sequence[ToricLatticeElement])-> ToricLatticeElement:
        return self.N([i/3  for i in (-self.K + sum(exceptional_curves))])
        # we can divide by 3 if we provide all exceptional curves for contracting to P2

    def _add_curve_name(self, curve:'Curve', name:str)->'Curve':
        curve.set_immutable()
        self._curve_names[curve] = name
        return curve

    def _curve_name(self, curve:'Curve') -> str:
        #print('naming curve ', curve)
        if curve in self._curve_names.keys() and curve in self.minus_one_curves+self.minus_two_curves:
            return self._curve_names[curve]
        else:
            return str(curve)

    def _curves_name(self, curves:Sequence['Curve']) -> str:
        return ', '.join([self._curve_name(c) for c in curves])#print('naming curve ', curve)

    @cached_property
    def minus_two_curves(self) -> list['Curve']:
        collinear = [self._add_curve_name(self.L - self.E[i] - self.E[j] - self.E[k], f'L_{i+1}{j+1}{k+1}') for i,j,k in self.collinear_triples]
        
        infinitesimal = [self._add_curve_name(self.E[chain[i]]-self.E[chain[i+1]], f'E_{chain[i]}{chain[i+1]}') for chain in self.infinitesimal_chains for i in range(len(chain)-1)]
        
        conic = [self._add_curve_name(2*self.L - sum(self.E[i] for i in six), f'Q{"".join(i+1 for i in range(len(self.E)) if i not in six)}') for six in self.sixs_on_conic]
        cubic = [3*self.L - sum(self.E) - self.E[i] for i in self.cusp_cubics]
        curves = collinear + infinitesimal + conic + cubic
        curves = normalize_rays(curves, self.N)
        return curves

    # def curve_name(self, c) -> str|None:
    #     if c in self.minus_one_curves:
    #         match self.dot(c, self.L):
    #             case 0:
    #                 pass



    @cached_property
    def minus_one_curves(self) -> list['Curve']:
        exceptional_curves = [self._add_curve_name(e, f'E_{i+1}') for i,e in enumerate(self.E)]
        lines = [self._add_curve_name(self.L-self.E[i]-self.E[j], f'L_{i+1}{j+1}') for i,j in itertools.combinations(range(len(self.E)), 2 )]
        conics = [self._add_curve_name(2 *self.L-sum(self.E[i] for i in points), f'Q_{"".join(i+1 for i in range(len(self.E)) if i not in points)}') for points in itertools.combinations(range(len(self.E)), 5 )]
        cubics = [3 *self.L-sum(points)-double for points in itertools.combinations(self.E, 7 ) for double in points]
        curves = exceptional_curves + lines + conics + cubics
        if self.degree == 1 :
            quartics = [4 *self.L-sum(self.E) - sum(double_points) for double_points in itertools.combinations(self.E, 3 )]
            quintics = [5 *self.L-sum(self.E) - sum(double_points) for double_points in itertools.combinations(self.E, 6 )]
            sextics = [6 *self.L-2 *sum(self.E) - triple_point for triple_point in self.E]
            curves += quartics + quintics + sextics
        #for c in curves:
         #   c.set_immutable()
        curves = normalize_rays(curves, self.N)
        if self.is_weak:
            curves = [c for c in curves if all(self.dot(c,f)>=0 for f in self.minus_two_curves)]
        return curves
        
    def curves_not_meeting(self, curves_to_filter, test_curves):
        return [c for c in curves_to_filter if all(self.dot(c,c2)==0  for c2 in test_curves)]

    @cached_property    
    def double_intersections(self): 
        '''
        Double intersections, which may coincide in triple ones (then 3 double intersections = 1 triple one)
        '''
        return [Point(frozenset([a,b])) for a,b in itertools.combinations(self.minus_one_curves, 2) if self.dot(a,b)>0 ]

    @cached_property    
    def triple_intersections(self): 
        '''
        Possible triple intersections
        '''
        return [Point(frozenset([a,b,c])) for a,b,c in itertools.combinations(self.minus_one_curves, 3 ) if self.dot(a,b)*self.dot(b,c)*self.dot(a,c)>0 ]

    #TODO refactor module to avoid usage of Ample in favor of NE through dualities. Reason is, NE has 240 rays in degree 1, and Ample has around 20k rays.
    @cached_property
    def Ample(self):
        return self.dual_cone(self.NE)

    # def independent_sets(self, curves, size = None):
    #     if size == None:
    #         yield from self.independent_sets(curves, 9 -self.degree)
    #         return
    #     if size == 0 :
    #         yield []
    #         return
    #     for i, v in enumerate(curves):
    #         orthogonals = [v2 for v2 in curves[i+1 :] if self.dot(v, v2)==0 ]
    #         for subset in self.independent_sets(orthogonals, size-1 ):
    #             yield subset + [v]

    def disjoint_subsets(self, curves:list, independent_with:list|None=None, maximal_only=True):
        '''
        
        EXAMPLES::
            sage: from collections import Counter
            sage: S = Surface(5)
            sage: Counter(len(s) for s in S.disjoint_subsets(S.minus_one_curves)).most_common()
            [(3, 10), (4, 5)]
        '''
        if independent_with == None:
            independent_with = []
        curves = [c for c in curves if all(self.dot(c,i)==0  for i in independent_with)]
        if len(curves) == 0 :
            yield []
            return
        curve = curves[-1]
        for subset in self.disjoint_subsets(curves[:-1], independent_with=independent_with):
            if all(self.dot(curve,c)==0  for c in subset):
                yield subset + [curve]
                if not maximal_only:
                    yield subset
            else:
                yield subset                

    def cone_representative(self, cone_type:str) -> 'NE_SubdivisionCone':
        return NE_SubdivisionCone.representative(self, cone_type)

    def cylinders(self, E, construction:str) -> Iterable['Cylinder']:
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
            case 'lines2':
                if self.degree <= 5:
                    for e1,e2,e3,e4 in itertools.combinations(E,int(4 )):
                        yield Cylinder.make_type_lines2(self, E, e1, e2, e3, e4)
                        yield Cylinder.make_type_lines2(self, E, e1, e3, e2, e4)
                        yield Cylinder.make_type_lines2(self, E, e1, e4, e2, e3)
            case 'tangent':
                if self.degree >= 3:  
                    for e in E:    
                        yield Cylinder.make_type_tangent(self, E, E_on_conic=[f for f in E if f!=e], E_on_tangent=[e])
            case 'P1xP1': 
                if self.degree <= 7  and self.degree >= 3 :
                    for Ei,Ej in itertools.permutations(E, int(2 )):
                        yield Cylinder.make_type_P1xP1(self, E, Ei, Ej)
            case 'cuspcubic':
                if self.degree == 2 :
                    for E_small in itertools.combinations(E,int(4 )):
                        yield Cylinder.make_type_cuspcubic(self, E, E_small)

    def all_cylinders(self, constructions):
        # TODO adapt the size of E for arbitrary degree
        for E in self.disjoint_subsets(self.minus_one_curves):
            for construction in constructions:
                yield from self.cylinders(E, construction)
    
    def __str__(self) -> str:
        return f"del Pezzo surface of degree {self.degree}"


class Contraction():    
    '''
    contracted_curves: A tuple of contracted (irreducible) curves.
    E: a list of exceptional reducible (-1)-curves, which form an orthogonal basis of the span of the contracted curves. They are sums of chains of contracted curves ending with a (-1)-curve.
    C: a list of contracted (irreducible) curves in the same order as E.
    '''
    S:Surface
    contracted_curves: tuple[ToricLatticeElement,...]
    E: tuple[ToricLatticeElement,...]
    C: tuple[ToricLatticeElement,...]

    def __init__(self, S:Surface, contracted_curves: Sequence[ToricLatticeElement]) -> None:
        """
        Initializes a new instance of the class.

        Args:
            S (Surface): The surface object.
            contracted_curves (Sequence[ToricLatticeElement]): A sequence of contracted curves.

        Returns:
            None
        """
        self.S = S
        self.contracted_curves = tuple(contracted_curves)
        E=[]
        C=[]
        assert self.is_valid(), 'contraction is not valid'
        for chain in self.chains_of_contracted_curves:
            for i in range(len(chain)):
                C.append(chain[i])
                E.append(sum(chain[i:]))
        assert len(E) == len(C)
        self.E = tuple(E)
        self.C = tuple(C)

    @cached_property
    def chains_of_contracted_curves(self) -> list[list[ToricLatticeElement]]:
        '''
        returns list of lists of curves of form [C1,C2..,Cn], where Ei are irreducible contracted curves such that Cn^2=-1, Ci^2=-2 for i<n, and Ci*C_{i+1}=1.
        '''
        minus_one_curves = [c for c in self.contracted_curves if self.S.dot(c,c)==-1]
        minus_two_curves = [c for c in self.contracted_curves if self.S.dot(c,c)==-2]
        def find_chain(minus_one_curve:ToricLatticeElement)->list[ToricLatticeElement]:
            chain = [minus_one_curve]
            while True:
                next_curves = [c for c in minus_two_curves if self.S.dot(c,chain[-1])==1 and c not in chain]
                match len(next_curves):
                    case 0:
                        return list(reversed(chain))
                    case 1:
                        chain.append(next_curves[0])
                    case _:
                        raise ValueError("more than one next curve, should be a chain")
        return [find_chain(c) for c in minus_one_curves]
                
    def _E_to_C(self, e:ToricLatticeElement)->ToricLatticeElement:
        return self.C[self.E.index(e)]

    def _C_to_E(self, c:ToricLatticeElement)->ToricLatticeElement:
        assert c in self.C, f"{c} is not in {self.C}"
        assert len(self.E) == len(self.C), f"{self.E} and {self.C} have different lengths"
        return self.E[self.C.index(c)]

    @cached_property
    def L(self):
        '''
        returns the class of line on the contracted surface (we assume it to be P2)
        '''
        assert self.is_P2()
        return self.S.Line(self.E)


    def project(self, v):
        return v + sum(self.S.dot(e,v)*e for e in self.E)
        #return self.S.N(self.picard_projection_matrix*v)

    def is_valid(self) -> bool:
        g = self.S.gram_matrix(self.contracted_curves)
        q = QuadraticForm(QQ,g)
        return q.is_negative_definite() and abs(g.det())==1

    def is_maximal(self) -> bool:
        return all(self.dot(c,c)>=0 for c in 
                   self.S.minus_one_curves+self.S.minus_two_curves)        

    def is_P2(self) -> bool:
        return self.is_maximal() and len(self.contracted_curves) == 9-self.S.degree

    def dot(self, a, b):
        return self.S.dot(self.project(a),self.project(b))

    def gram_matrix(self, rays):
        return matrix([[self.dot(a,b) for a in rays] for b in rays])

    def __str__(self) -> str:
        return f"contraction of curves {self.contracted_curves} on {self.S}"

    def make_cylinder_QT(self, E_on_conic:Sequence[ToricLatticeElement], E_on_tangent:Sequence[ToricLatticeElement], E_on_fibers: list[Sequence[ToricLatticeElement]]|None = None)->'Cylinder':  
        return Cylinder.make_type_tangent_weak(self.S, self, E_on_conic, E_on_tangent, E_on_fibers)

    def find_cylinders_LL(self)->Generator['Cylinder',None,None]:
        for L1, dim_1 in self.line_classes:
            for L2, dim_2 in self.line_classes:
                intersection = [self._E_to_C(e) for e in self.E if self.S.dot(e,L1)==1 and self.S.dot(e,L2)==1]
                if len(intersection)>1:
                    continue
                base_curves = [chain[0] for chain in self.chains_of_contracted_curves]
                curves_on_other_fibers = [c for c in base_curves if self.S.dot(self._C_to_E(c),L1+L2)==0]
                complement = list(self.C)
                # we don't check the case when, say, fibers L12, L34 and L56 intersect in the same point (i.e., Echkardt point). Indeed, this case is automatically covered when generic flexibility is in question
                support = list(self.C) + [L1, L2] + [self.L- self._C_to_E(c) for c in curves_on_other_fibers]
                if dim_1 == 0:
                    complement.append(L1)
                if dim_2 == 0:
                    complement.append(L2)
                dimension = dim_1 + dim_2 if len(intersection) == 0 else 0
                transversal = dimension>0 

                basepoint = Point(frozenset([L1,L2])) if len(intersection)==0 else None
                yield Cylinder.make(self.S, 
                                         construction='lines',
                                         contraction=self,
                                         pencil_generators=(L1,L2),
                                         complement=complement, 
                                         support=support, 
                                         fiber=2*self.L - sum(intersection), 
                                         basepoint=basepoint,
                                         transversal=transversal,
                                         dimension=dimension)

    def find_cylinders_QT(self)->Generator['Cylinder',None,None]:
        labels = ['Q','T','QT','F1','F2','F3','F4']
        #TODO remove configs?
        empty_config = {c:'' for c in self.C}
        for Q, dim_Q in self.conic_classes:
            for T, dim_T in self.line_classes:
                config = empty_config.copy()
                QT_is_linear = True
                for c in self.C: 
                    e = self._C_to_E(c)
                    if self.S.dot(e,Q) == 1:
                        config[c] = 'Q'
                    if self.S.dot(e,T) == 1:
                        config[c] += 'T'
                intersection = [c for c,v in config.items() if v=='QT']                
                dim_QT = dim_Q + dim_T
                # if(dim_QT<=2):
                #     ic(Q,T, config, dim_Q, dim_T, intersection)
                if len(intersection)>2:
                    continue
                elif len(intersection)==2:
                    c1,c2 = intersection
                    if self.S.dot(c1,c2) != 1:
                        continue
                elif len(intersection) == 1:
                    c = intersection[0]
                    c_neighours = [d for d in self.C if self.S.dot(d,c)==1]
                    if len(c_neighours)==1 and config[c_neighours[0]]!='':
                            continue
                    dim_QT -= 1
                else:
                    dim_QT -= 1
                    QT_is_linear = False
                
                if dim_QT < 0:
                    continue

                transversal = dim_QT > 0 or not QT_is_linear
                complement = list(self.C)
                support = complement + [Q, T]                
                if dim_Q == 0:
                    complement.append(Q)
                    if dim_QT == 0 and QT_is_linear:
                        complement.append(T)
                elif dim_T == 0:
                    complement.append(T)
                #TODO special fibers

                basepoint = Point(frozenset([Q,T])) if len(intersection)==0 else None
                yield Cylinder.make(self.S, 
                                         complement=complement, 
                                         support=support, 
                                         construction='tangent',
                                         contraction=self,
                                         pencil_generators=(Q,T),
                                         fiber=2*self.L - sum(intersection), 
                                         basepoint=basepoint,
                                         transversal=transversal,
                                         dimension=dim_QT)

    def are_collinear(self, c1, c2, c3) -> bool:
        E = [self._C_to_E(c) for c in (c1,c2,c3)]
        return self.L-sum(E) in self.S.minus_two_curves

    @cached_property
    def conic_classes(self) -> list[tuple[ToricLatticeElement,int]]:
        '''
        return a list of Picard classes of proper transforms of conics in P2, each with the dimension of the corresponding linear system.
        '''
        result = []
        for p in itertools.product(*[range(len(chain)+1) for chain in self.chains_of_contracted_curves]):
            dimension = 5 - sum(p)
            if dimension < 0:
                continue
            curves_on_conic = [curve for i,n in enumerate(p) for
                               curve in self.chains_of_contracted_curves[i][:n]]
            if any(self.are_collinear(c1,c2,c3) for c1,c2,c3 in itertools.combinations(curves_on_conic,3)):
                continue
            conic_class = 2*self.L - sum(self._C_to_E(c) for c in curves_on_conic)
            conic_class.set_immutable()
            result.append([conic_class, dimension])
        return result
    
    @cached_property
    def line_classes(self) -> list[tuple[ToricLatticeElement,int]]:
        '''
        return a list of Picard classes of proper transforms of lines in P2, each with the dimension of the corresponding linear system.
        UPD: actually, they are linear systems, i.e., smaller objects than classes
        '''
        #TODO rewrite by taking all collinear triples, then non-collinear pairs, then single base curves, and then a generic class
        result = []
        for p in itertools.product(*[range(len(chain)+1) for chain in self.chains_of_contracted_curves]):
            dimension = 2 - sum(p)
            if sum(p)>3:
                continue
            curves_on_line = [c for i,n in enumerate(p) for c in self.chains_of_contracted_curves[i][:n]]
            if any(self.are_collinear(c1,c2,c3)==False for c1,c2,c3 in itertools.combinations(curves_on_line,3)):
                continue
            dimension = 2 - len(curves_on_line)
            if dimension == -1:
                if self.are_collinear(*curves_on_line):
                    dimension = 0
                else:
                    continue
            
            # check that if two ci are on line, then another ci collinear with them is on line
            other_curves_meeting_line = [c for c in self.C if any(self.are_collinear(c1, c2, c) for c1, c2 in itertools.combinations(curves_on_line, 2))]            
            if not set(other_curves_meeting_line).issubset(set(curves_on_line)):
                continue
            
            line_class = self.L - sum(self._C_to_E(c) for c in curves_on_line)
            line_class.set_immutable()
            result.append([line_class, dimension])
        return result



    #def check_restraints_on_line(self, curves:Sequence[ToricLatticeElement]) -> bool:




class Curve(ToricLatticeElement):
    #pass #TODO find a reasonable implementation
    # TODO make immutable    
    def __init__(self, coordinates, name:str|None=None):
        super().__init__(coordinates)
        self.name = name
    #components = []

    def __mul__(self, other:'Curve') -> int:
        product = self.S._N_to_M(other) * self
        return product

    


@dataclass(frozen=True)
class Point:
    curves : frozenset[Curve]
    # TODO define a proper equality

@dataclass(frozen=True)
class Cylinder:
    '''
    fiber is the generic fiber class of the cylinder (or cylinders if several). If the collection is transversal, then the fiber is empty.
    complement is the list of (-1)-curves lying in the complement of the cylinders' union
    support is the list of (-1)-curves lying in the support of divisors, where the cylinders' union is polar
    '''
    S : Surface
    construction: str|None
    contraction: Contraction|None
    pencil_generators: tuple[Curve, Curve]|None
    complement : tuple[Curve, ...]
    support : tuple[Curve, ...]
    fiber: Curve|None = None
    basepoint: Point|None = None
    transversal: bool|None = None
    dimension: int|None = None

    @classmethod
    def make(cls, S:Surface, complement:list[Curve], support:list[Curve], fiber:Curve|None=None, basepoint:Point|None=None, construction:str|None=None, contraction:Contraction|None=None, pencil_generators:tuple[Curve, Curve]|None=None, transversal:bool|None=None, dimension:int|None=None) -> 'Cylinder':
        for c in complement:
            c.set_immutable()
        for c in support:
            c.set_immutable()
        if fiber!=None:
            fiber.set_immutable()    
        return cls(S, construction, contraction, pencil_generators, tuple(complement), tuple(support), fiber, basepoint, transversal, dimension)

    @classmethod
    def make_type_lines(cls, S:Surface, E:list[Curve], e:Curve)->'Cylinder':
        '''
        We draw lines through image of e and image of one of E, for each of E.
        '''
        L = S.Line(E)
        complement = E+[L-e-f for f in E if f!=e]
        return cls.make(S, complement, complement, L-e, construction='lines', transversal=False)

    @classmethod
    def make_type_lines2(cls, S:Surface, E:list[Curve], e1:Curve, e2:Curve, e3:Curve, e4:Curve)->'Cylinder':
        '''
        See Perepechko-2013, Section 3.1.
        e1, e2, e3, e4 are distinct members of E.
        We draw lines L1 through images of e1 and e2 in P2, L2 through images of e3 and e4. Their images provide a pencil of lines, which gives a cylinder. Other members of E may be assumed to be in different singular fibers.
        '''
        others = [e for e in E if e!=e1 and e!=e2 and e!=e3 and e!=e4]
        assert len(others) == len(E)-4 
        L = S.Line(E)
        lines = [L-e1-e2, L-e3-e4]
        for l in lines:
            l.set_immutable()
        complement = E + [L-e1-e2, L-e3-e4] + [L-e for e in others]
        return cls.make(S, complement, complement, L, basepoint=Point(frozenset(lines)), construction='lines2', transversal = False)

    @classmethod
    def make_type_tangent(cls, S:Surface, E:list[Curve], E_on_conic:Sequence[Curve], E_on_tangent:Sequence[Curve], E_on_fibers: Sequence[Sequence[Curve]]|None = None)->'Cylinder':
        '''
        We draw a conic and a tangent on P2, and specify the positions of blown up points.

        
        Special case:
        See Cheltsov-Park-Won Ex.4.1.2 and Perepechko Ex.5.2. 
        We draw a conic through all blown up points but one, and a tangent line through the remaining one. 
        Technically, for each e we obtain a collection of two cylinders for two tangents respectively. 
        This pair has no common fibers.
        '''
        if E_on_fibers == None:
            E_on_fibers = []
        disjoint_subsets_of_E = [set(E_on_conic), set(E_on_tangent)] + [set(E_on_fiber) for E_on_fiber in E_on_fibers]
        assert len(E) == sum(len(e) for e in disjoint_subsets_of_E) 
        assert len(E) == len(set.union(*disjoint_subsets_of_E))
        L = S.Line(E)
        tangent = L - sum(e for e in E_on_tangent)
        conic = 2*L - sum(e for e in E_on_conic)
        support = E + [conic, tangent] + [2*L - sum(E_on_fiber) for E_on_fiber in E_on_fibers]

        # checks if there are at least two tangents through a given point. 
        two_tangents_through_point = len(E_on_tangent)<=1 and all(len(E_on_fiber)<=1 for E_on_fiber in E_on_fibers)<=1

        #fibers containing at least two Ei
        special_fibers = [E_on_fiber for E_on_fiber in E_on_fibers if len(E_on_fiber)>1]

        # there are at least two tangents such that there is a fiber containing two given points, see CPW Th.6.2.3
        two_tangents_with_fiber_through_two_points = len(E_on_tangent)==0 and len(special_fibers)<=1 and all(len(E_on_fiber)<=2 for E_on_fiber in special_fibers)

        if two_tangents_through_point or two_tangents_with_fiber_through_two_points:            
            complement = E + [conic]
            transversal = True
        else:
            complement = support
            transversal = None
        return cls.make(S, complement, support, 2*L, construction='tangent', transversal=transversal)
    
    @classmethod
    def make_type_tangent_weak(cls, S:Surface, contraction:Contraction, E_on_conic:Sequence[Curve], E_on_tangent:Sequence[Curve], E_on_fibers: list[Sequence[Curve]]|None = None)->'Cylinder':
        '''
        We draw a conic and a tangent on P2, and specify the positions of blown up points.

        Parameters:
            S: the surface
            contraction: the chosen contraction to P2, i.e. a list of contracted curves.
            E_on_Q: the minus-curves on conic Q
            E_on_T: the minus-curves on tangent T. May intersect with E_on_conic.
            E_on_fibers: for each fiber containing at least two minus-curves (except for the intersection of Q and T), there is a list of minus-curves on that fiber. Fibers containing one minus-curve can be omitted. The minus-curves on the intersection of Q and T are omitted.
        
        Remark:
            Instead of (-2)-curves we operate with the corresponding reducible (-1)-curve, i..e., the total transform of a (-1)-curve. So, if E2 is inf. near of E1, and conic preimage contains E2, then it contains E1.
        '''
        if E_on_fibers == None:
            E_on_fibers = []
        
        # TODO just provide classes of L, Q, T, Fi?
        #TODO assertions: lines are in T and not in Q
        assert all(e in contraction.E for e in itertools.chain(E_on_conic,E_on_tangent,*E_on_fibers))

        L = S.Line(contraction.E)
        tangent = L - sum(e for e in E_on_tangent)
        conic = 2*L - sum(e for e in E_on_conic)
        fibers = [2*L - sum(E_on_fiber) for E_on_fiber in E_on_fibers]
        #fibers containing at least two Ei
        special_fibers = [E_on_fiber for E_on_fiber in E_on_fibers if len(E_on_fiber)>1]

        for e in contraction.E:
            if S.dot(e,tangent)==0 and S.dot(e,conic)==0 and all(S.dot(e,fiber)==0 for fiber in E_on_fibers):
                fibers.append(2*L - e)
        pol_gens = list(contraction.contracted_curves) + [conic, tangent] + fibers
        assert all(S.dot(c,e)<=1 for c in pol_gens for e in contraction.E)


        # the dimension of the linear system of possible conics Q
        dim_Q = 5 - len(E_on_conic)
        # the dimension of the linear system of possible lines T (before the tangency condition)
        dim_T = 2 - len(E_on_tangent)

        if dim_Q<0 or dim_T<0:
            raise ValueError('No cylinders exist')

        # now we compute parameters for the moduli of (Q,T) with tangency condition

        intersection_QT = [e for e in E_on_conic if e in E_on_tangent]
        if len(intersection_QT)>2:
            raise ValueError('Q and T have more than double intersection point')
        elif len(intersection_QT)==2:
            QT_is_linear = True
            QT_dimension = dim_Q + dim_T
        else:
            QT_is_linear = False
            QT_dimension = dim_Q + dim_T - 1

        for fiber in special_fibers:
            QT_dimension -= len(fiber) - 1
            QT_is_linear = False

        if QT_dimension<0:
            raise ValueError('No cylinders exist')     
        transversal = QT_dimension > 0 or not QT_is_linear
        forb_gens = pol_gens if transversal else [c for c in pol_gens if c!=tangent]

        # checks if there are at least two tangents through a given point. 
        two_tangents_through_point = len(E_on_tangent)<=1 and all(len(E_on_fiber)<=1 for E_on_fiber in E_on_fibers)<=1
        # there are at least two tangents such that there is a fiber containing two given points, see CPW Th.6.2.3
        two_tangents_with_fiber_through_two_points = len(E_on_tangent)==0 and len(special_fibers)<=1 and all(len(E_on_fiber)<=2 for E_on_fiber in special_fibers)
        if two_tangents_through_point or two_tangents_with_fiber_through_two_points:            
            assert transversal == True

        return cls.make(S, forb_gens, pol_gens, 2*L, construction='tangent', transversal=transversal, dimension = QT_dimension)

    @classmethod
    def make_type_tangent2(cls, S:Surface, E:list[Curve], E_on_conic:Sequence[Curve])->'Cylinder':
        '''
        See Perepechko-2013 Section 4.1. 
        We draw a conic through blown up points at E_small and take a generic tangent line. 
        Technically, we obtain a collection of cylinders with no common fibers.
        '''
        assert len(E_on_conic)<=5
        return cls.make_type_tangent(S, E, E_on_conic, [], [[e] for e in E if e not in E_on_conic])


    @classmethod
    def make_type_P1xP1(cls, S:Surface, E:list[Curve], Ei,Ej)->'Cylinder':
        '''
        Cheltsov-Park-Won Example 4.1.6 and Lemma 4.2.2 for contraction of E1..E4, L-E5-E6 in degree 3. See also Perepechko Ex. 5.3.
        This is again a pair of cylinders without common fibers.
        #TODO can we use it for degree 2?
        '''
        assert S.degree<=7  and S.degree>=3 
        assert Ei in E and Ej in E and Ei!=Ej
        L = S.Line(E)
        E_complement = [e for e in E if (e!=Ei) and (e!=Ej)]
        conic = -S.K-L+Ei
        complement = [e for e in E_complement]+[L-Ei-Ej, conic]
        support = complement + [L-Ei,L-Ej]
        return cls.make(S, complement, support, None, construction='P1xP1')

    @classmethod
    def make_type_cuspcubic(cls, S:Surface, E:list[Curve], E_small:Sequence[Curve])->'Cylinder':
        '''
        See Cheltsov-Park-Won Ex.4.1.13 and Th.6.2.2 case 2.
        C is an anticanonical cuspidal curve, cubic on contraction to P^2 defined by (E minus E_small) and M.
        We assume that there are at least two such curves (with distinct cusp points).
        E_small defines the contraction from degree 5.
        '''
        assert S.degree >= 2  and S.degree <= 5
        assert len(E_small) == 4 and all(e in E for e in E_small)
        L = S.Line(E)
        C = - S.K
        M = [2*L-sum(E_small)] + [L-e for e in E_small]        
        complement =  [e for e in E if e not in E_small] 
        support = [C] + M + complement
        fiber = 2*C
        return cls.make(S, complement, support, fiber, construction='cuspcubic', transversal=True)


    @classmethod
    def make_type_cuspcubic2(cls, S:Surface, E:list[Curve], E_small:Sequence[Curve])->'Cylinder':
        #TODO adapt for degree <=5
        '''
        See Cheltsov-Park-Won Ex.4.1.11.
        '''
        assert S.degree == 2 
        assert len(E_small) == 4 
        L = S.Line(E)
        C = 3 *L-sum(E)
        M = [2 *L-sum(E_small)] + [L-e for e in E_small]
        support = [C] + M + E
        complement = E
        return cls.make(S, complement, support, None, construction='cuspcubic')

    @classmethod
    def make_type_cuspcubic_L12(cls, S:Surface, E:list[Curve], E_on_line:Sequence[Curve])->'Cylinder':
        '''
        See Cheltsov-Park-Won Th.6.2.2, Case 2.4
        '''
        assert S.degree >= 2 and S.degree <= 7
        assert len(E_on_line) == 2 and all(e in E for e in E_on_line)
        E_other = [e for e in E if e not in E_on_line]
        L = S.Line(E)
        C = 3 *L-sum(E)
        L12 = L - sum(E_on_line)
        F = [L-e for e in E_on_line]
        complement = E_other + [L12]
        support = complement + [C] + F
        return cls.make(S, complement, support, 2*(F[0]+F[1]), construction='cuspcubic')


    @cached_property
    def Pol(self):
        return Cone(self.support)

    @cached_property
    def Forb(self):
        return Cone(self.complement)

    def is_polar_on(self, cone:ConvexRationalPolyhedralCone|str):
        if isinstance(cone, str):
            cone = NE_SubdivisionCone.representative(self.S, cone)
        return relint_contains_relint(self.Pol, cone)
    
    def is_complete_on(self, cone:ConvexRationalPolyhedralCone, exclude:ConvexRationalPolyhedralCone|None=None):
        '''
        checks if the collection is H-complete for ample divisor classes H from the relative interior of cone
        exclude is a cone of divisors to be excluded from completeness check
        '''
        intersection = cone.intersection(self.Forb)
        if not relint_contains_relint(cone, intersection):
            return True
        if exclude == None:
            return False
        forb_intersection_excluded = all(exclude.contains(ray) for ray in intersection.rays())
        return forb_intersection_excluded

    def is_transversal(self):
        if self.transversal!=None: 
            return self.transversal 
        return self.fiber == None

    def compatible_representatives(self, complete=False) -> Iterable[str]:
        '''
        yields types of representatives, on which self is polar and, optionally, complete
        '''
        types = NE_SubdivisionCone.cone_types(self.S)
        for t in types:
            #print(f'looking at type {t}')
            cone = NE_SubdivisionCone.representative(self.S, t)
            if self.is_polar_on(cone):
                if complete and not self.is_complete_on(cone):
                    continue
                yield t

    def latex(self) -> str:
        result = \
            f'''
            construction: {self.construction}
            used contraction: {self.S._curves_name(self.contraction.contracted_curves)}
            pencil generators: {self.S._curves_name(self.pencil_generators)}
            generic fiber: {self.fiber}
            polarity cone generators: {self.S._curves_name(self.support)}
            forbidden cone generators: {self.S._curves_name(self.complement)}
            is transversal: {self.transversal}
            '''
            
        return result


class CylinderList(list):
    '''
    A list of cylinders with extra info.

    complement is the list of (-1)-curves lying in the complement of the cylinders' union
    fiber is the generic fiber class of the collection. If the collection is transversal, then the fiber is empty
    '''
    def __init__(self, cylinders:Iterable[Cylinder], S:Surface|None=None) -> None:
        super().__init__(cylinders)
        if len(self)==0  and S==None:
            raise ValueError('The surface is not defined')
        self.S = S if S!=None else self[0].S
        for cyl in self:
            self._validate_cylinder(cyl)

    def __add__(self, other):        
        if self.S != other.S:
            raise ValueError('Different surfaces')
        cylinders = list(self) + list(other)
        return CylinderList(cylinders, self.S)
    
    def __setitem__(self, index, item):
        super().__setitem__(index, self._validate_cylinder(item))

    def insert(self, index, item):
        print(super(), type(self))
        super().insert(index, self._validate_cylinder(item))

    def append(self, item):
        super().append(self._validate_cylinder(item))

    def extend(self, other):
        if isinstance(other, type(self)):
            if self.S != other.S:
                raise ValueError('Different surfaces')
            super().extend(other)
        else:
            super().extend(self._validate_cylinder(item) for item in other)

    def _validate_cylinder(self, cylinder) -> Cylinder:
        if not isinstance(cylinder, Cylinder):
            #print(type(cylinder), cylinder.S)
            raise TypeError(f'{cylinder} is not a Cylinder')
        if cylinder.S != self.S:
            raise ValueError(f'{cylinder} is defined on other surface')
        return cylinder

    @cached_property
    def fiber(self):
        fiber = self[0].fiber
        if fiber!=None:
            if any(fiber!=cyl.fiber for cyl in self[1:]):
                return None
        return fiber
        
        
    @cached_property
    def complement(self):
        if len(self)==0:
            raise ValueError('The collection is empty, complement curves are undefined')
            return self.S.minus_one_curves
        complement = [curve for curve in self[0].complement if all(curve in cyl.complement for cyl in self[1:])]
        return complement
        
    @cached_property
    def complement_double_points(self):
        return [p for p in self.S.double_intersections if all(p in curve for curve in self.complement)]
    
    @cached_property
    def Pol(self):
        #TODO ensure that the dimension of cone does not change
        result = Cone([],lattice=self.S.N.dual()).dual() # ambient space
        for cylinder in self:
            result = result.intersection(cylinder.Pol)
        interior_element = sum(result.rays())
        if all(c.Pol.relative_interior_contains(interior_element) for c in self):
            return result
        else:
            return Cone([],lattice=self.S.N)

    @cached_property #TODO if it is rarely invoked, then we should not cache the result. Are there tests for memory usage?
    def Forb(self):
        if len(self) == 0 :
            return Cone([],lattice=self.S.N.dual()).dual()
        return Cone(self.complement, lattice=self.S.N)

    def is_polar_on(self, cone):
        '''
        cone is either a Cone or a type of cone, which would be converted to a representative
        '''
        if isinstance(cone, str):
            cone = NE_SubdivisionCone.representative(self.S, cone)
        return relint_contains_relint(self.Pol, cone)

    def make_polar_on(self, cone) -> 'CylinderList':
        '''
        cone is either a Cone or a type of cone, which would be converted to a representative
        '''
        if isinstance(cone, str):
            cone = NE_SubdivisionCone.representative(self.S, cone)
        cylinders = [c for c in self if c.is_polar_on(cone)]
        return CylinderList(cylinders, self.S)

    def reduce(self, keep_double_points:bool=False) -> 'CylinderList':
        '''
        keep_double_points is a boolean that indicates whether double points in the union of cylinders should be preserved
        '''
        if keep_double_points:
            raise NotImplementedError
        cylinders_count = Counter()
        removed_cylinders = set()
        for cylinder in self:
            if all(cylinders_count[curve]>0  for curve in cylinder.complement):
                removed_cylinders.add(cylinder)
                continue
            cylinders_count.update(cylinder.complement)
        while True:
            for cyl in self:
                if cyl in removed_cylinders:
                    continue
                if all(cylinders_count[curve]>1  for curve in cyl.complement):
                    removed_cylinders.add(cyl)
                    for curve in cyl.complement:
                        cylinders_count[curve]-=1 
                    break
            else:
                break
        return CylinderList([c for c in self if c not in removed_cylinders], self.S)
        #TODO implement this for double (triple) points, some old code below
        if len(complement_points(collection)['double'])!=0 :
            print('this collection does not cover all double points')
            return collection
        print('deleting cylinders, remaining:')
        while True:
            for i in range(len(collection)):
                collection_without_cylinder = collection[:i] + collection[i+1 :]
                if len(complement_points(collection_without_cylinder)['double'])==0 :
                    del collection[i]
                    print(len(collection), end=' ')
                    break
            else:
                print('no cylinders to delete')
                break
        return collection

    def is_transversal(self):
        if any(cyl.transversal == True for cyl in self):
            return True
        if self.fiber == None:
            return True
        if len(self) == 0:
            return False
        basepoint = self[0].basepoint
        if basepoint == None:
            return False
        return any(basepoint != cyl.basepoint for cyl in self[1:])

    
    def is_complete_on(self, cone:ConvexRationalPolyhedralCone, exclude:ConvexRationalPolyhedralCone|None=None):
        '''
        checks if the collection is H-complete for ample divisor classes H from the relative interior of cone
        exclude is a cone of divisors to be excluded from completeness check
        '''
        intersection = cone.intersection(self.Forb)
        if not relint_contains_relint(cone, intersection):
            return True
        if exclude == None:
            return False
        forb_intersection_excluded = all(exclude.contains(ray) for ray in intersection.rays())
        return forb_intersection_excluded


    def is_generically_flexible_on(self, cone, exclude=None, restrict_to_ample=True):
        '''
        checks if the collection provides generic flexibility for ample divisors in the provided cone
        exclude is a cone of divisors to be excluded from completeness check
        '''
        if isinstance(cone, str):
            cone = NE_SubdivisionCone.representative(self.S, cone)
        if isinstance(exclude, str):
            exclude = NE_SubdivisionCone.representative(self.S, cone)
        if restrict_to_ample:
            cone = cone.intersection(self.S.Ample)
            if not relint_contains_relint(self.S.Ample, cone):
                return True
        is_polar = self.is_polar_on(cone)
        is_complete = self.is_complete_on(cone, exclude)
        return is_polar and is_complete and self.is_transversal()
    
    def __str__(self) -> str:
        return f'Cylinder collection with {len(self)} cylinders on {self.S}'
        

class NE_SubdivisionCone(ConvexRationalPolyhedralCone):
    '''
    This class represents a cone obtained by subdivisions from NE and keeps the creation info.
    '''

    @classmethod
    def NE(cls, S:Surface):
        if S.is_weak:
            curves = S.minus_one_curves + S.minus_two_curves
        else:
            curves = S.minus_one_curves
        NE = cls(curves, S.N)
        NE.S = S
        NE.parent = None
        NE.parent_face = None
        NE.parent_ray = None
        return NE
        

    @classmethod
    def from_face_and_ray(cls, parent:'NE_SubdivisionCone', face:ConvexRationalPolyhedralCone, ray:ToricLatticeElement):
        assert face.is_face_of(parent)
        assert face != parent
        rays = list(face.rays())+[ray]
        rays = sorted(normalize_rays(rays, parent.S.N))
        child = cls(rays, lattice=parent.S.N)
        child.S = parent.S 
        child.parent = parent
        child.parent_face = face
        child.parent_ray = ray
        return child
    
    
    def subdivision_faces(self, ray: ToricLatticeElement) -> Iterable[ConvexRationalPolyhedralCone]:
        '''
        returns faces of self that do not contain ray
        '''
        for faces_of_dim in self.faces()[:-1]:
            for face in faces_of_dim:       # type: ignore
                if not face.contains(ray):  # type: ignore
                    yield face              # type: ignore
        
    def _subdivision_ray(self) -> ToricLatticeElement:
        '''
        returns the standard dividing ray of the cone depending on its type
        '''
        match self.type:
            case 'NE':
                return -self.S.K
            case 'C':
                return sum(self.parent_face.rays())
            case _:
                return sum(self.rays())


    def children(self, subdivision_ray: ToricLatticeElement | None = None)->Iterable['NE_SubdivisionCone']:
        '''
        returns a generator of subdivision cones
        '''
        if subdivision_ray == None:
            subdivision_ray = self._subdivision_ray()
        for f in self.subdivision_faces(subdivision_ray):
            cone = NE_SubdivisionCone.from_face_and_ray(self, f, subdivision_ray)
            if not relint_contains_relint(self, cone):
                continue
            yield cone


    def ample_children(self, subdivision_ray: ToricLatticeElement | None = None)->Iterable['NE_SubdivisionCone']:
        '''
        returns a generator of subdivision cones whose relative interiors contain ample divisors
        
        EXAMPLES ::
        sage: cone_C = next(c for c in Surface(5).NE.ample_children() if c.type == 'C')
        sage: sorted(Counter(c.type for c in cone_C.ample_children()).items())
        [('C(0)', 1), ('C(1)', 6), ('C(2)', 12), ('C(3)', 4), ('C(3)-P1xP1', 4)]
        '''
        if subdivision_ray == None:
            subdivision_ray = self._subdivision_ray()
        ray_is_ample = self.S.Ample.relative_interior_contains(subdivision_ray)
        for cone in self.children(subdivision_ray):
            if ray_is_ample or cone.has_ample_divisors():
                yield cone
            


    def make_child(self, face:ConvexRationalPolyhedralCone, subdivision_ray: ToricLatticeElement | None = None)->'NE_SubdivisionCone':
        '''
        returns a subdivision cone corresponding to the given face
        '''
        # TODO convert ValueErrors into return None
        if subdivision_ray == None:
            subdivision_ray = self._subdivision_ray()
        if not self.contains(subdivision_ray):
            raise ValueError(f'make_child: subdivision_ray {subdivision_ray} is not contained in this cone ({self} with rays {list(self.rays())}) of type {self.type}')
        if not face.is_face_of(self):
            raise ValueError(f'make_child: face (rays {list(face.rays())}) is not a face of this cone ({self} with rays {list(self.rays())}) of type {self.type}')
        cone = NE_SubdivisionCone.from_face_and_ray(self, face, subdivision_ray)
        if not relint_contains_relint(self, cone):
            raise ValueError(f'make_child: resulted cone (rays {list(cone.rays())}, type {cone.type}) does not intersect the relative interior of this cone ({self} with rays {list(self.rays())}) of type {self.type}')
        return cone
        # should we check for containing ample divisors? maybe not 

    def print(self):
        info = f'Cone of type {self.type}, rays\n{self.rays()}, on {self.S}'
        if self.parent != None:
            info +=f' with parent of type {self.parent.type}, parent ray {self.parent_ray}, parent face rays {list(self.parent.rays())}.'
        return info
        
    



    @cached_property
    def Ample(self):
        '''
        returns a cone, which relative interior consists of all ample classes in self.cone
        '''
        #TODO we assume that L is in NE, make this explicit
        ample = self.intersection(self.S.Ample)
        if self.S.Ample.relative_interior_contains(sum(ample.rays())):
            return ample
        else:
            return Cone([],lattice=self.S.N)

    @property 
    def Ample2(self): # is this faster than Ample?
        dual_rays = self.S.NE.rays() + self.S.dual_cone(self).rays()
        return self.S.dual_cone(Cone(dual_rays))

    def has_ample_divisors(self):
        return self.Ample.dimension() > 0
        
    @cached_property
    def type(self):
        '''
        Returns the Fujita type of the parent face for children of NE.
        EXAMPLES::
            sage: from collections import Counter
            sage: Counter(c.type for c in Surface(4).NE.ample_children()).most_common()
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
        elif self.parent == self.S.NE and self.parent_ray == -self.S.K:
            rays = list(self.parent_face.rays())
            if not all(r in self.S.minus_one_curves for r in rays):
                return 'NE.unknown-child-1'
            E = next(self.S.disjoint_subsets(rays))
            if len(E) == len(rays):
                res = f'B({rank})'
                other_disjoint_curves = self.S.curves_not_meeting(self.S.minus_one_curves, rays)
                if rank == 9 - self.S.degree - 1 and len(other_disjoint_curves)==0:
                    res+='-P1xP1'
                return res
            # below is not a complete check for type C. Do we need one? yes!
            meeting_curves = [[r for  r in rays if self.S.dot(r,e)==1 ] for e in E]
            if not all(len(curves)==1  for curves in meeting_curves):
                return 'NE.unknown-child-2'
            E1 = [curves[0] for curves in meeting_curves]
            B = E[0]+E1[0]
            if all(E[i]+E1[i]==B for i in range(len(E))) and \
                all(self.S.dot(e1,e2)<=0 for e1 in E1 for e2 in E1) and len(E)==rank-1:
                return 'C'
            return 'NE.unknown-child-3'
        elif self.parent.type == 'C' and self.parent_ray == self.parent._subdivision_ray():
            rays = list(self.parent_face.rays())
            B = self.parent_ray
            E = [r for r in rays if r in self.S.minus_one_curves]
            has_anticanonical_ray = -self.S.K in rays
            if len(rays) - len(E) != 1 if has_anticanonical_ray else 0:
                #print(rays, E, B, rank)
                return 'C.unknown-child-6'
            if not all(self.S.dot(e1,e2)<=0 for e1 in E for e2 in E):
                return 'C.unknown-child-4'
            if not all(self.S.dot(B,e)==0 for e in E):
                return 'C.unknown-child-5'
            cone_type = f'C({len(E)})'
            if len(E) == len(self.S.E) - 1:
                if all(any(self.S.dot(c,e)!=0 for e in E) for c in self.S.minus_one_curves):
                    cone_type += '-P1xP1'
            if not has_anticanonical_ray:
                cone_type += '-no-K'
            return cone_type
        return self.parent.type + '.unknown-child-0'

    #this method works slowly in degree 1, why?
    @classmethod 
    @cache 
    def representative(cls, S:Surface, cone_type:str):
        parent = S.NE
        match list(cone_type):
            case ['N','E']:
                return S.NE
            case ('B','(', num_rays, ')'):
                num_rays = int(num_rays)
                assert num_rays == int(re.search(r'\d+', cone_type).group()) # type: ignore
                face_rays = S.E[:num_rays]
            case ['C']:
                face_rays = S.E[:-1 ] + [S.L - e - S.E[-1 ] for e in S.E[:-1 ]]
            case ('C','(', num_rays, ')'):
                num_rays = int(num_rays)
                assert num_rays < 9 - S.degree
                face_rays = S.E[:num_rays] + [-S.K]
                parent = NE_SubdivisionCone.representative(S, 'C')
            case ('C','(', num_rays, ')', '-','n','o','-','K'):
                num_rays = int(num_rays)
                assert num_rays < 9 - S.degree
                face_rays = S.E[:num_rays]
                parent = NE_SubdivisionCone.representative(S, 'C')
            case _:
                num_rays = len(S.E) - 1
                if cone_type == f'B({num_rays})-P1xP1':
                    face_rays = S.E[:-2]+[S.L - S.E[-1] - S.E[-2]]
                elif cone_type == f'C({num_rays})-P1xP1':
                    face_rays = S.E[:-2] + [S.L - S.E[-1] - S.E[-2]] + [-S.K]
                    parent = NE_SubdivisionCone.representative(S, 'C')
                elif cone_type == f'C({num_rays})-P1xP1-no-K':
                    face_rays = S.E[:-2] + [S.L - S.E[-1] - S.E[-2]]
                    parent = NE_SubdivisionCone.representative(S, 'C')
                else:
                    raise ValueError('unknown cone type', cone_type)
        face_rays = sorted(normalize_rays(face_rays, S.N))
        face = Cone(face_rays, lattice=S.N)
        return parent.make_child(face)

    @staticmethod
    def cone_types(S:Surface):
        blowups = 9 - S.degree
        return ['NE', f'B({blowups-1 })-P1xP1', f'C({blowups-1 })-P1xP1', 'C']+            [f'B({i})' for i in range(blowups+1 )] + [f'C({i})' for i in range(blowups)]
        #removed:  [f'C({blowups-1 })-P1xP1-no-K']+ [f'C({i})-no-K' for i in range(blowups)]

    @classmethod
    def representatives(cls,S:Surface):
        for cone_type in cls.cone_types(S):
            yield cls.representative(S, cone_type)

    def relative_volume(self, cone:ConvexRationalPolyhedralCone):
        '''
        given a subcone of self, computes the intersection with an affine hyperplane orthogonal to -K and returns the volume of the subcone divided by the volume of self
        '''
        # see https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/base7.html#sage.geometry.polyhedron.base7.Polyhedron_base7.volume
        raise NotImplementedError

