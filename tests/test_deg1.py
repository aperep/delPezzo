from delPezzo import Surface, CylinderList, Cylinder, NE_SubdivisionCone
import pytest

'''
We follow [CPW15]_
'''

@pytest.mark.skip(reason="needs a rewrite")
def test_deg1():
    S1 = Surface(1)
    E = S1.E
    L = S1.L
    C1 = Cylinder.make_type_tangent(S1, E, E[3:], [E[2]], [[e] for e in E[:2]])
    pol_C1 = E + [2*L-E[0]-E[1], L-E[2], 2*L-sum(E[3:])]
    for r in pol_C1:
        r.set_immutable()
    assert set(pol_C1) == set(C1.Pol.rays())
    assert CylinderList([C1]).is_generically_flexible_on(C1.Pol)
    return
    C1_compatible_types = C1.compatible_representatives(complete=True)
    C1_compatible_types_expected = ['B(3)', 'C(3)']

    assert set(C1_compatible_types_expected) == set(C1_compatible_types), f'compatible representatives should be {C1_compatible_types_expected}, got {C1_compatible_types}'

    C2 = Cylinder.make_type_cuspcubic_L12(S2, S2.E, S2.E[-2:])
    C2_compatible_types_expected = ['B(6)-P1xP1', 'C(6)-P1xP1']
    C2_compatible_types = C2.compatible_representatives(complete=True)
    assert set(C2_compatible_types) == set(C2_compatible_types_expected), f'compatible representatives should be {C2_compatible_types_expected}, got {C2_compatible_types}'
