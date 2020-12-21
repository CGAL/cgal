* Use the exact tag in the Delaunay triangulation.
* After converting exact to inexact, we may have equal points in the bbox faces. We should fix that. Can we?
* Can we accelerate initializer when using exact kernel?
* In the merge pvertices event, check for parallel cases not only for the cropping case but for the other cases, too!
* I should not use pface_of_pvertex() but rather find the correct face that will intersect or not the corresponding iedge and if I should decrement k for it then I do, otherwise I stop/continue. I am pretty sure that pface_of_pvertex() is wrong especially for the open case, back || front subcase.


1. The first step: we compute N full line arrangements. That gives an intersection graph. But to get it, we need to intersect 3D planes. Right?
2. How do they handle cases with multiple coplanar polygons?
3. Can we avoid kinetic completely? What if we use labeling instead?
4. What is the exact stop criteria for k intersections? When two polygons intersect at the very beginning, does it count as an intersection? How to count it when we meet an iedge or an ivertex? Can we squeeze through the whole or not?
5. It looks like for open case, we never insert new polygons. Or better, do we need open case at all? Maybe it should be reassigned into closing case! Because, in the original paper, they do not have such a case.
6. Graph-cut based surface extraction does not guarantee a 2-manifold. What if we simply stop when we meet a roof polygon or another wall polygon? In this case, we do not guarantee convex polyhedra but who cares?
7. Do we use the trick from the paper when inserting only the first three events (corresponding to the three closest lines) of an active vertex in the queue instead of all possible events?
8. What about multiple simultaneous collisions between k primitives at the same time? Do we handle that? Do we need to add a random perturbation as in the original code?
9. How to do the subdivision and merging efficiently?
10. How should we choose the pface whose k is updated? Can it differ for different event types?
11. Do we need to add extra triangles as in my experimental code to handle k?


QUESTIONS:
1. Do we need to intersect all 3D planes in order to find N full line arrangements?
- Yes.

2. How do you handle cases with multiple coplanar polygons?
- Merge polygons, which are at the same cell.

3. Can we avoid kinetic completely?
- Probably yes, but there is a problem of how to find the correct criteria that guarantees the volume convexity. But we can use the trick with tagging pfaces at least for all initial polygons.

4. Do you have any guarantee on the number of final volumes? Is this number optimal for each k?
- No, it is super difficult.

5. When two polygons intersect at the very beginning, does it count as an intersection? Can we squeeze through the whole or not?
- Both are ok.

Ideas:
Do the functions collision_occured/is_occupied work correctly? Not sure.
Should I use uniform k with respect to the support plane rather than pface? I think yes.
Try to use for pvertex->iedge events the function version: is_occupied(pvertex, ivertex, iedge).
Study again the case with 7 and 12 input polygons. They have 1 and 3 hanging pfaces respectively for k = 1. 27 k = 1 for building b and 19 k = 2 for the same building.
Try to simplify the case as much as possible such that choice of pface does not change anything, we do not have many events until the bad case happens, we are sure that we correctly passed all intermediate events, use only the case k = 1 so that we do not propagate.
Try to do again what I already tried in the experimental code but for the new refactored case.
Try to use strict version of the function is_occupied() everywhere.
What if I extend the initially intersected polygons until they fill the whole cell, then hanging pfaces for the case of the test with 7 and 12 polygons (test-8 and test-9) should not happen!
Try to see what JP is doing in more detail and compare his steps with ours!
Where I should use which version of is_occupied()?
