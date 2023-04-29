import sympy as sp
import itertools
# lines on a cubic https://www.geogebra.org/geometry/ynm2c5y9

column_sum = lambda M: sum((M.col(i) for i in range(M.cols)), sp.Matrix([0]*M.rows))

K = sp.Matrix([-3, 1, 1, 1, 1, 1, 1]) # canonical divisor
Line = lambda exceptionals : (-K + column_sum(exceptionals)) / 3

Q = sp.diag(1, -1, -1, -1, -1, -1, -1) # intersection form in basis L, E1, .., E6
def meet(D1, D2):
  return (D2.T*Q*D1)[0,0]
#incidence_matrix = Matrix([[meet(i,j) for j in minus_one_curves] for i in minus_one_curves])


def minus_one_curves():
    '''
    returns a list of all 27 (-1)-curves on the cubic surface
    '''
    E = sp.eye(7)[:,1:]
    exceptional = [E.col(i) for i in range(E.cols)] # exceptional curves are basis vectors e2...e7    
    L = Line(E)

    lines = [L-ei-ej for ei,ej in itertools.combinations(exceptional, 2)]
    conics = [-K-L+e for e in exceptional]
    return exceptional + lines + conics


def independent_sets(vectors, size = 6):
  if size == 0:
    yield []
  for i, v in enumerate(vectors):
    orthogonals = [v2 for v2 in vectors[i+1:] if meet(v, v2)==0]
    for subset in independent_sets(orthogonals, size-1):
      yield subset + [v]


    
def cylinder_cone_rays(E):
  '''
  E is a sympy matrix with 6 exceptional curves as columns
  returns rays of the cone of divisors H such that a cylinder (P^2 minus cubic C and tangent line) is H-polar
  '''
  tangent = Line(E) - E[:, -1]
  conic = - K - tangent
  return sp.Matrix.hstack(E, tangent, conic)

  

def cylinders():
    '''
    returns a list of cones of compatible divisors for cylinders of considered type
    '''
    cylinders = []
    for blowdown in independent_sets(minus_one_curves()):
      line = Line(sp.Matrix.hstack(*blowdown))
      for e in blowdown:
        cylinders.append(blowdown+[line-e,-K-line+e])
    return [[tuple(d) for d in c] for c in cylinders]



