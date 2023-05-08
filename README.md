# delPezzo

This is a Sagemath module for studying flexibility of affine cones over del Pezzo surfaces.

To use it, install Sagemath and put the file `delPezzo_cylinders.py` in your working directory. The following example code imports the module,
creates a Surface instance, chooses the subdivision cone B3 and a cylinder collection.
```
from delPezzo_cylinders import *
S = Surface(3)
B3 = S.cone_representative('B(3)')
collection = Cylinder.make_type_cuspcubic(S, S.E, S.E[-4:])
```

The list of available subdivision cones for `S` is returned by `NE_SubdivisionCone.cone_types(S)`.
The rays of polarity and forbidden cones of the collection are returned by `collection.Pol.rays()`
and `collection.Forb.rays()` respectively.
The properties of the collection in the relative interior of a given cone (B3 here) as well as subdivision cones, 
where the collection is polar and complete, can be checked with the following methods.
```
collection.is_polar_on(B3) # False
collection.is_complete_on(B3) # True
collection.is_transversal() # True
collection.is_generically_flexible_on(B3) # False
list(collection.compatible_representatives()) # ['B(2)', 'C(2)']
list(collection.compatible_representatives(complete=True)) # ['B(2)', 'C(2)']
```

The method `S.all_cylinders(constructions)` returns a collection comprised of
all cylinders of certain constructions (e.g., `constructions = ['lines','tangent']`)
corresponding to all choices of the contraction to the projective plane. Such a collection is useful in conjuction with the
following methods.

The method `collection.make_polar_on(cone)` filters out the cylinders that are
not polar inside the cone, and `collection.reduce()` removes abundant cylinders from
the collection while keeping the forbidden cone unchanged.
