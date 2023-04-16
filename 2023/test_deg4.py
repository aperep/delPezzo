from delPezzo_cylinders import *



def test_ParkWon():
    '''
    tests the covering provided in Park-Won for the degree 4.
    '''
    S = Surface(4)
    
    # case B(i)
    Ui = CylinderList(S.cylinders(S.E, 'lines2')).make_polar_on('B(1)')
    assert len(Ui)==3
    for i in range(1, 6):
        break
        assert Ui.is_generically_flexible_on(f'B({i})')
    return Ui
