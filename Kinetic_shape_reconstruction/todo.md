* When we create a new face, should we use the initial k for it, e.g. 2 or should we copy the last counted k, e.g. 2 or 1?
* Should we decrease k in the polygon splitter when two polygons are intersected at the very beginning?
* Should we count the number of intersections per polygon or per mesh face?
* Should we keep facei-facej number of intersections?
* Use the exact tag in the Delaunay triangulation.
* After converting exact to inexact, we may have equal points in the bbox faces. We should fix that.
* Can we accelerate initializer when using exact kernel?