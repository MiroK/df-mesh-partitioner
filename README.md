# Partitioning of FEniCS meshes by PyMetis

In a nutshell
```python
import dolfin as df
from partitioning import *

mesh = df.UnitSquareMesh(32, 32)
graph = mesh_2_nxgraph(mesh, edge_weight_f=lambda u, v: u+v)

cell_f = partition(graph, 5, weighted=False)
df.File('results/unweighted.pvd') << cell_f

cell_f = partition(graph, 5, weighted=True)
df.File('results/weighted.pvd') << cell_f
```

## Dependencies
- FEniCS 2019.1.+  (python3)	
- [`pymetis`](https://github.com/inducer/pymetis/blob/main/pymetis/__init__.py)

## TODO
- [ ] Do I need node weights?