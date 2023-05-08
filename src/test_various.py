from delPezzo_cylinders import *

def test_S5_covering():
    S = Surface(5)
    constructions = ['P1xP1', 'lines', 'tangent', 'lines2']
    cylinders = [c for c in S.all_cylinders(constructions)]
    total_collection = CylinderList(cylinders)
    cones = [c for c in NE_SubdivisionCone.representatives(S)]
    cone = NE_SubdivisionCone.representative(S,'B(0)')
    collection = CylinderList(cylinders)
    covering = collection.make_polar_on(cones[1 ]).reduce()
    assert len(covering) == 5 
    return covering

def test_cone_types(deg:int=3):
    S = Surface(deg)
    types = NE_SubdivisionCone.cone_types(S)
    for t in types:
        assert t == NE_SubdivisionCone.representative(S,t).type, f'types {t} and {NE_SubdivisionCone.representative(S,t).type} are not same'


def test_S3_tangent_covering():
    S = Surface(3)
    collection = CylinderList(S.cylinders(S.E, 'tangent'))
    for i in range(1, 7):
        t = f'B({i})'
        cone = NE_SubdivisionCone.representative(S, t)
        assert len(collection.make_polar_on(cone).reduce())>0 
    cone = NE_SubdivisionCone.representative(S, 'B(0)')
    assert len(collection.make_polar_on(cone).reduce())==0 

def test_S3_cover_anticanonical():
    S = Surface(3)
    constructions = ['P1xP1', 'lines', 'tangent']
    cylinders = [c for c in S.all_cylinders(constructions)]
    cone = NE_SubdivisionCone.representative(S,'B(0)')
    collection = CylinderList(cylinders)
    covering = collection.make_polar_on(cone)
    assert len(covering) == 0


if __name__=="__main__":
    S = Surface(3)

