from delPezzo_cylinders import Surface


def test_blowup_P2():
    P2 = Surface(9)
    blowups = list(P2.blowups())
    assert len(blowups) == 1
    assert blowups[0].degree == 8


def test_blowup_deg_8():
    S = Surface(8)
    blowups = list(S.blowups())
    assert len(blowups) == 2


def test_contraction_P2():
    S = Surface(8)
    blowups = list(S.blowups())
    b,c = blowups
    if b.is_weak:
        b,c = c,b
    assert not b.is_weak
    contractions = list(c.contractions_P2())
    assert len(contractions) == 1
    contraction = contractions[0]
    assert len(contraction.contracted_curves) == 2
    assert b.N([0,0,1]) in contraction.contracted_curves
    assert contraction.is_valid()
    assert contraction.is_maximal()
    assert all(contraction.project(c)==contraction.S.N(0) for c in contraction.contracted_curves)