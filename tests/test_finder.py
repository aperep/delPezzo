from delPezzo import Finder, Surface


def test_iterate():
    S = Surface(5)
    F = Finder(S)
    assert F.iterate()
    assert len(F.coverings) == 1
    assert len(F.coverings[0].collection) == 3