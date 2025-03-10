from delPezzo import Surface, CylinderList, Cylinder, NE_SubdivisionCone, WeakDependencies


'''
We follow Perepechko 2023
'''

S1 = Surface(6, )

def test_deg6_weak_1():
    S = Surface(6, dependencies=WeakDependencies(collinear_triples=[[1,2,0]]))
    cylinder = Cylinder.make(S, 
                    complement=S.E,
                    support=S.E + [S.L-e for e in S.E],
                    fiber=S.L,
                    transversal=True
                    )
    assert cylinder.is_complete_on(S.Ample)
    assert cylinder.is_polar_on(S.Ample)
    assert cylinder.is_transversal()


def test_deg6_weak_2():
    S = Surface(6, dependencies=WeakDependencies(infinitesimal_chains=[[0,1]]))
    E = S.E
    L = S.L
    # a-e are exceptional curves that comprise a zigzag
    a,b,c,d,e = L-E[0]-E[1], E[1], E[0]-E[1], L-E[0]-E[2], E[2]    
    cyl1 = Cylinder.make(S, 
                    complement=[b,c,d,e],
                    support=[b,c,d,e]+[L],
                    fiber=S.L,
                    transversal=True
                    )
    Lprime = 2*L-sum(E)
    cyl2 = Cylinder.make(S, 
                complement=[a,b,c,d],
                support=[a,b,c,d]+[Lprime],
                fiber=Lprime,
                transversal=True
                )
    col = CylinderList([cyl1, cyl2])
    assert col.is_generically_flexible_on(S.Ample)

