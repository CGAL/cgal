* Use the exact tag in the Delaunay triangulation.
* After converting exact to inexact, we may have equal points in the bbox faces. We should fix that. Can we?
* Can we accelerate initializer when using exact kernel?
* In the merge pvertices event, check for parallel cases not only for the cropping case but for the other cases, too!
* I should not use pface_of_pvertex() but rather find the correct face that will intersect or not the corresponding iedge and if I should decrement k for it then I do, otherwise I stop/continue. I am pretty sure that pface_of_pvertex() is wrong especially for the open case, back || front subcase.