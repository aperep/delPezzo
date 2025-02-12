import re
from typing import Generator
from delPezzo import Surface
from importlib import resources

def generate_surfaces(degree:int) -> Generator[Surface, None, None]:
    results = []
    pattern_list = r"\[\d+(?:[ \t]*,[ \t]*\d+)+\]"
    with resources.files("delPezzo").joinpath("Lubbes_list.txt").open() as file: 
        lines = file.readlines()
        for line in lines:
            index, rank, line = line.strip().split(' ',2)
            if rank != str(9-degree):
                continue
            list_contents, root_type = line.split(']',1)
            basis_labels = re.findall(r'\d+', list_contents)
            collinear_triples = [[int(i)-1 for i in b[1:]] for b in basis_labels if 
                                len(b)==4 and b[0]=='1']
            
            infinitesimal_chains = []
            for i,j in [b for b in basis_labels if len(b)==2]:
                i,j = int(i)-1, int(j)-1
                if len(infinitesimal_chains) == 0:
                    infinitesimal_chains.append([i,j])
                else:
                    if infinitesimal_chains[-1][-1] == i:
                        infinitesimal_chains[-1].append(j)
                    else:
                        infinitesimal_chains.append([i,j])
            
            sixs_on_conic = [[i for i in range(8) if i+1 not in (int(b[1]),int(b[2]))] for b in basis_labels if len(b)==3 and b[0]=='2']
            cusp_cubics = [int(b[2])+1 for b in basis_labels if len(b)==3 and b[:2]=='30']
            extra = {'Lubbes_type':root_type, 'Lubbes_basis':basis_labels, 'Lubbes_index':index}
            yield Surface(degree, collinear_triples=collinear_triples, infinitesimal_chains=infinitesimal_chains, sixs_on_conic=sixs_on_conic, cusp_cubics=cusp_cubics, extra=extra)
                

if __name__ == "__main__":
    for surface in generate_surfaces(5):
        print(surface)