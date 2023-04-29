from delPezzo_cylinders import Surface, CylinderList, Cylinder, NE_SubdivisionCone



def test_Perepechko2020():
    '''
    tests the coverings provided in Perepechko, 2020 for the degree 3.
    '''
    S = Surface(3)
    
    U1 = CylinderList([Cylinder.make_type_tangent(S, S.E, E_on_conic=S.E[1:] , E_on_tangent=[S.E[0]])])
    U2 = CylinderList([Cylinder.make_type_tangent(S, S.E, E_on_conic=S.E[:-1], E_on_tangent=[S.E[-1]])])
    def V(i,j):
        return Cylinder.make_type_P1xP1(S, S.E, S.E[i-1], S.E[j-1])
    U3 = CylinderList([V(5,6)])
    U4 = CylinderList([V(5,6), V(6,5)])
    U5 = CylinderList([V(5,6), V(4,6)])

    def representatives_where_collection_is_flexible(U):
        return set(t for t in NE_SubdivisionCone.cone_types(U.S) if U.is_generically_flexible_on(t))

    assert representatives_where_collection_is_flexible(U1) == set(f'B({i})' for i in range(1,7))
    assert representatives_where_collection_is_flexible(U2) == {'B(5)', 'B(6)', 'C(5)'}
    assert representatives_where_collection_is_flexible(U4) == {'B(5)-P1xP1', 'C(5)-P1xP1'}
    assert representatives_where_collection_is_flexible(U5) == {f'C({i})' for i in range(0,4)}

    C4 = NE_SubdivisionCone.representative(S, 'C(4)')
    forb_dual_ray = S.L + S.E[-2]-2*S.E[-1]
    C4prime = S.dual_cone([forb_dual_ray,-forb_dual_ray]).intersection(C4)
    assert len(C4prime.rays()) == 5
    assert U3.is_generically_flexible_on(C4, exclude=C4prime)
    assert U5.is_generically_flexible_on(C4prime)

