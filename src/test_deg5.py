from delPezzo_cylinders import Surface, CylinderList



def test_Perepechko2013():
    '''
    tests the covering provided in Perepechko, 2013 for the degree 5.
    '''
    S = Surface(5)
    cylinders = CylinderList(list(S.all_cylinders(['lines2'])))
    assert cylinders.is_generically_flexible_on(S.Ample)
    return cylinders
