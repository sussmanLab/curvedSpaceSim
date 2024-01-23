# Change log

### Changes in progress

* submesher: switch to breadth-first search with early stopping (optional, not obviously better)

* submesher: check aggresive submeshing routines for potential improvements

### version 0.9

* basic interface between MD-like framework and CGAL's Xin and Wang single-source-all-points geodesic algorithm implemented
* cell-list spatial sorting and Submeshing routines to greatly improve performance in the case of short-ranged forces implemented
* MPI functionality implemented
