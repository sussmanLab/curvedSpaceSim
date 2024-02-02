# Citations for curvedSpaceSimulation

If you use this software packagefor a publication or project, please cite the main curvedSpaceSimulation methods paper:

(1) (arXiv and eventual journal link here)

The core algorithm used to compute geodesic distances and paths comes from:

## Xin and Wang exact discrete geodesic distance (implemented in triangleMeshSpace)

(2) [Shi-Qing Xin and Guo-Jin Wang](https://dl.acm.org/doi/10.1145/1559755.1559761). Improving chen and han's algorithm on the discrete geodesic problem. ACM Trans. Graph., 28(4):104:1–104:8, September 2009

This paper is itself an improvement on the  work of 

(3) Chen and Han (Shortest paths on a polyhedron. Internat. J. Comput. Geom. Appl., 6:127–144, 1996),

which itself improves on the original MMP algorithm for finding discrete geodesics:

(4) Joseph S. B. Mitchell, D. M. Mount, and C. H. Papadimitriou. The discrete geodesic problem. SIAM J. Comput., 16:647–668, 1987

The Xin and Wang algorithm is as-implemented by CGAL:

(5a) CGAL,Computational Geometry Algorithms Library, http://www.cgal.org
(5b) Stephen Kiazyk and Sébastien Loriot and Éric Colin de Verdière. Triangulated Surface Mesh Shortest Paths. In CGAL User and Reference Manual. CGAL Editorial Board, 5.6 edition, 2023

## Vector heat method (implemented in goemetryCentralMeshSpace)

(2) [Nicholas Sharp, Yousuf Soliman, and Keenan Crane](https://dl.acm.org/doi/abs/10.1145/3243651). The Vector Heat Method. ACM Trans. Graph 38(3):1-19 (2019)

This is built into the Geometry-central framework, which has the following bibtex info:
@article{geometrycentral,
  title={GeometryCentral: A modern C++ library of data structures and algorithms for geometry processing},
  author={Nicholas Sharp and Keenan Crane and others},
  howpublished="\url{https://geometry-central.net/}",
  year={2019}
}
