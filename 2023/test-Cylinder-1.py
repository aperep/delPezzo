
# Generated by CodiumAI

import pytest

"""
Code Analysis

Main functionalities:
The Cylinder class represents a cylinder on a surface, which is a collection of curves that do not intersect and whose union is a divisor on the surface. The class provides methods to compute the polar and forbidden cones of the cylinder, and to check if a given cone is polar with respect to the cylinder.

Methods:
- __init__: initializes a Cylinder object with a Surface object, a list of curves representing the complement of the cylinder, and an optional curve representing the generic fiber of the cylinder.
- is_polar_on: checks if a given cone is polar with respect to the cylinder.
- Pol: returns the polar cone of the cylinder.
- Forb: returns the forbidden cone of the cylinder.

Fields:
- S: a Surface object representing the surface on which the cylinder is defined.
- complement_curves: a list of Curve objects representing the curves lying in the complement of the cylinder.
- fiber: an optional Curve object representing the generic fiber of the cylinder.
"""



class TestCylinder:

    # Tests that a cylinder object can be initialized with a surface object, a list of curve objects representing the complement of the cylinder, and an optional curve object representing the generic fiber of the cylinder. tags: [happy path]
    def test_cylinder_creation(self):
        S = Surface(_sage_const_2)
        C = Cylinder(S, [S.E[_sage_const_0], S.E[_sage_const_1]], S.L)
        assert C.S == S
        assert C.complement_curves == [S.E[_sage_const_0], S.E[_sage_const_1]]
        assert C.fiber == S.L

    # Tests that the forb method computes the forbidden cone of the cylinder correctly. tags: [general behavior]
    def test_forbidden_cone(self):
        S = Surface(_sage_const_3)
        C = Cylinder(S, [S.E[_sage_const_0], S.E[_sage_const_1], S.E[_sage_const_2]], S.L)
        assert C.Forb == Cone([S.N([-_sage_const_1, _sage_const_1, _sage_const_1]), S.N([_sage_const_1, -_sage_const_1, _sage_const_1]), S.N([_sage_const_1, _sage_const_1, -_sage_const_1])])

    # Tests that the pol method computes the polar cone of the cylinder correctly. tags: [general behavior]
    def test_polar_cone(self):
        S = Surface(_sage_const_4)
        C = Cylinder(S, [S.E[_sage_const_0], S.E[_sage_const_1], S.E[_sage_const_2]], S.L)
        assert C.Pol == Cone([S.N([-_sage_const_1, _sage_const_0, _sage_const_0]), S.N([_sage_const_0, -_sage_const_1, _sage_const_0]), S.N([_sage_const_0, _sage_const_0, -_sage_const_1]), S.N([_sage_const_1, _sage_const_1, _sage_const_1])])

    # Tests that a cylinder object can be initialized with an empty list of complement curves. tags: [edge case]
    def test_empty_complement_curves(self):
        S = Surface(_sage_const_4)
        c = Cylinder(S, [])
        assert c.complement_curves == []

    # Tests that a cylinder object cannot be initialized with a list of complement curves that do not lie in the complement of any cylinder on the surface. tags: [edge case]
    def test_invalid_complement_curves(self):
        S = Surface(_sage_const_4)
        with pytest.raises(ValueError):
            c = Cylinder(S, [S.L])

    # Tests that the is_polar_on method correctly checks if a given cone is polar with respect to the cylinder. tags: [general behavior]
    def test_is_polar_on(self):
        S = Surface(_sage_const_4)
        E = S.minus_one_curves
        c = Cylinder(S, E[:_sage_const_3], E[_sage_const_3])
        cone = Cone([E[_sage_const_3], E[_sage_const_4]])
        assert c.is_polar_on(cone) == True
        cone2 = Cone([E[_sage_const_3], E[_sage_const_5]])
        assert c.is_polar_on(cone2) == False