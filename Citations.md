# Citations for curvedSpaceSimulation

If you use this software packagefor a publication or project, please cite the main curvedSpaceSimulation methods paper:

(1) [The current arXiv version of the paper](https://arxiv.org/abs/2404.19751): Toler H. Webb and Daniel M. Sussman. "curvedSpaceSim: A framework for simulating particles interacting along geodesics." arXiv preprint arXiv:2404.19751 (2024).


The core algorithm used to compute geodesic distances and paths comes from

(2) [Shi-Qing Xin and Guo-Jin Wang](https://dl.acm.org/doi/10.1145/1559755.1559761). Improving chen and han's algorithm on the discrete geodesic problem. ACM Trans. Graph., 28(4):104:1–104:8, September 2009

This paper is itself an improvement on the  work of 

(3) Chen and Han (Shortest paths on a polyhedron. Internat. J. Comput. Geom. Appl., 6:127–144, 1996),

which itself improves on the original MMP algorithm for finding discrete geodesics:

(4) Joseph S. B. Mitchell, D. M. Mount, and C. H. Papadimitriou. The discrete geodesic problem. SIAM J. Comput., 16:647–668, 1987

The Xin and Wang algorithm is as-implemented by CGAL:

(5a) CGAL,Computational Geometry Algorithms Library, http://www.cgal.org
(5b) Stephen Kiazyk and Sébastien Loriot and Éric Colin de Verdière. Triangulated Surface Mesh Shortest Paths. In CGAL User and Reference Manual. CGAL Editorial Board, 5.6 edition, 2023
