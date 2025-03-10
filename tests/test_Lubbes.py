from delPezzo import generate_surfaces, generate_surface, CylinderGenerator

def test_generate_surfaces():
    list4 = list(generate_surfaces(4))
    list7 = list(generate_surfaces(7))
    assert len(list4) == 16
    assert len(list7) == 2

def test_generate_surface():
     line = "18 5 [23, 45] 2A1"
     S = generate_surface(line)
     assert len(S.minus_two_curves) == 2

def test_Surface_Line_Lubbes18():
     line = "18 5 [23, 45] 2A1"
     S = generate_surface(line)
     list(CylinderGenerator.all_cylinders(S, ['lines']))