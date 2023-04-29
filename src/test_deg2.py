from delPezzo_cylinders import Surface, CylinderList, Cylinder, NE_SubdivisionCone


'''
We follow CPW
'''

S2 = Surface(2)

def test_deg2_cuspidal():
    C1 = Cylinder.make_type_cuspcubic(S2, S2.E, S2.E[3:])
    C1_compatible_types = C1.compatible_representatives(complete=True)
    C1_compatible_types_expected = [f'B({i})' for i in range(3,8)] + [f'C({i})' for i in range(3,7)]

    assert set(C1_compatible_types_expected) == set(C1_compatible_types), f'compatible representatives should be {C1_compatible_types_expected}, got {C1_compatible_types}'

    C2 = Cylinder.make_type_cuspcubic_L12(S2, S2.E, S2.E[-2:])
    C2_compatible_types_expected = ['B(6)-P1xP1', 'C(6)-P1xP1']
    C2_compatible_types = C2.compatible_representatives(complete=True)
    assert set(C2_compatible_types) == set(C2_compatible_types_expected), f'compatible representatives should be {C2_compatible_types_expected}, got {C2_compatible_types}'
