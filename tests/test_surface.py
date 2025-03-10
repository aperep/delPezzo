from delPezzo import Surface
from collections import Counter
            
def test_disjoint_subsets():
    S = Surface(5)
    C = Counter(len(s) for s in S.disjoint_subsets(S.minus_one_curves)).most_common()
    assert C == [(3, 30), (2, 30), (1, 10), (4, 5), (0, 1)]