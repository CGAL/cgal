* Use the exact tag in the Delaunay triangulation.
* After converting exact to inexact, we may have equal points in the bbox faces. We should fix that. Can we?
* Can we accelerate initializer when using exact kernel?
* In the merge pvertices event, check for parallel cases not only for the cropping case but for the other cases, too!