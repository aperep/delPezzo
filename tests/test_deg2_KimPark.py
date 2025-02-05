from dataclasses import dataclass
from typing import Any
from delPezzo import Surface, CylinderList, Cylinder, NE_SubdivisionCone

S2 = Surface(2)


'''
this file is abandoned.

(move to txt/)

The notation of Kim-Park is as follows:

The anticanonical cuspidal curve is denoted by C and its singular point by p.
The anticanonical tacnodal pair is denoted by C1, C1', and their intersection point by q.

For the blowdown to P1xP1, we denote by L1 and L2 the (1,0) and (0,1) lines through p/q.

For the blowdown to F1=Bl(P2), we denote by E the (-1)-curve and by M the 0-line through p/q.

For the blowdown from dP6 to P2, we denote by T_0 the tangent to C at p, by T_1 the conic at P2 tangent to T_0 at p and passing through images of Ei. So, T0,T1 are 1-curves in S.

For the blowdown to dP5, we denote by Gi the 0-curves not passing through some of Ej (which ones?)

'''

@dataclass
class Kim_Park_Cylinder:
    number: int
    cone_types: list[str]
    cylinder: Cylinder
    surface_info: dict[str,Any]

    @classmethod
    def make(cls, number:int):
        match number:
            case 5:
                cone_types = ['B(3)']
                cylinder = Cylinder.make_type_cuspcubic(S2, S2.E, S2.E[3:])
                surface_info = {'type':'II'}
            case _:
                raise NotImplementedError
                cone_types = ['B()']
                cylinder = Cylinder.make_type_()
                surface_info = {'type':'II'}
        return cls(number, cone_types, cylinder, surface_info)



def test_KimPark_II():
    assert CylinderList([Kim_Park_Cylinder.make(5).cylinder]).is_generically_flexible_on('B(3)')
    pass


        

