# delPezzo

This is a Sagemath module for studying flexibility of affine cones over del Pezzo surfaces, see [arXiv:2305.06462](https://arxiv.org/abs/2305.06462). 



## Quickstart

The module `delPezzo_cylinders.py` is intended to be used with a Sagemath-provided Python interpreter `sage -python`.

1. install Sagemath locally or in a container, e.g., using a VSCode [Dev Container](https://code.visualstudio.com/docs/devcontainers/containers) extension with the included `.devcontainer/devcontainer.json` file.
2. Put the file `src/delPezzo_cylinders.py` in your working directory.
3. Import the module with `from delPezzo_cylinders import *` in your script or a Jupyter notebook.
4. Run your code in one of the following ways:
   - Python script `script.py` with `sage -python script.py`
   - Sagemath script `script.sage` with `sage script.sage`
   - Jupyter notebook with a kernel provided by Sagemath, see [instructions](https://doc.sagemath.org/html/en/installation/launching.html), or in case of the Dev Container, select Python interpreter `sage -python`.

## Usage

The following example code imports the module,
creates a Surface instance, chooses the subdivision cone **B3** and a cylinder collection.

```python
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

```python
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


## Tests

The tests are implemented in files `src/test_*.py` and can be mass-checked with `pytest`. There is also a notebook `tests/tests.ipynb` with all tests included for convenience.