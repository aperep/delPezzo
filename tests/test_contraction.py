from delPezzo import Surface, Contraction, generate_surface

def test_chains_of_contracted_curves():
    S = generate_surface("18 5 [23, 45] 2A1")
    contraction = S.standard_contraction()
    curve1, curve2 = contraction.chains_of_contracted_curves[1]
    assert S.dot(curve1, curve2) == 1
    assert S.dot(curve1, curve1) == -2
    assert S.dot(curve2, curve2) == -1
    pass

def test_valid():
    S = generate_surface("18 5 [23, 45] 2A1")
    contraction = S.standard_contraction()
    assert contraction.is_valid()
    assert Contraction.make_valid_contraction(S, S.E) != None
