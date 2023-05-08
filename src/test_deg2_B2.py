from sage.geometry.cone import Cone 
from sage.geometry.polyhedron.constructor import Polyhedron
import tikzplotlib

from delPezzo_cylinders import Surface, Cylinder, NE_SubdivisionCone, cones_intersect_by_relint


S2 = Surface(2)

#cone = NE_SubdivisionCone.representative(S2, 'B(2)')
#cone = Cone([-S2.K+S2.E[0]+S2.E[1]])

E_on_conic = S2.E[2:]
E1, E2 = S2.E[0], S2.E[1]
C1 = Cylinder.make_type_tangent(S2, S2.E, E_on_conic, [E1], [[E2]])
C2 = Cylinder.make_type_tangent(S2, S2.E, E_on_conic, [E2], [[E1]])
C3 = Cylinder.make_type_tangent(S2, S2.E, E_on_conic, [], [[E1, E2]])

#B2_cone = NE_SubdivisionCone.representative(S2, 'B(2)').intersection(S2.Ample)
#B2_cone_section = Polyhedron(vertices=[-S2.K, -S2.K + S2.E[0], -S2.K + S2.E[0] + S2.E[1], -S2.K + S2.E[1]])

def point_to_square(point):
    x,y = point[1:3]
    return x+1, y+1

def polyhedron_to_square(P):
    vertices = [point_to_square(v) for v in P.vertices()]
    return Polyhedron(vertices=vertices)

def square_point_to_ray(S:Surface,x,y):
    x1, x0 = x.numerator(), x.denominator()
    y1, y0 = y.numerator(), y.denominator()
    x, y, z = x0*y0, x1*y0, x0*y1
    return -x*S.K + y*S.E[0] + z*S.E[1]
    

def polar_subset_on_B2(S:Surface, cone):
    vertices=[-S.K, -S.K + S.E[0], -S.K + S.E[0] + S.E[1], -S.K + S.E[1]]
    B2_cone = Cone(vertices)
    B2_cone_section = Polyhedron(vertices=vertices)
    P = cone.polyhedron().intersection(B2_cone_section)
    P = polyhedron_to_square(P)
    if P.is_empty():
        return Polyhedron(vertices=[(0,0)])
    if cones_intersect_by_relint(cone, B2_cone):
        return P
    return Polyhedron(vertices=[(0,0)])



if __name__ == '__main__':
    #P0 = polar_subset_on_B2(S2, B2_cone)
    P1 = polar_subset_on_B2(S2, C1.Pol)
    P2 = polar_subset_on_B2(S2, C2.Pol)
    P3 = polar_subset_on_B2(S2, C3.Pol)
    #P0.plot() + P1.plot() + P2.plot() + P3.plot()
    P = P1.plot() + P2.plot() + P3.plot()
    with open('../local/B2_CPW.tex', 'w') as f:
        print(tikzplotlib.get_tikz_code(P.matplotlib()), file=f)
