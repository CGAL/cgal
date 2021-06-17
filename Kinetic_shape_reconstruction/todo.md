* Use the exact tag in the Delaunay triangulation.
* After converting exact to inexact, we may have equal points in the bbox faces. We should fix that. Can we?
* Can we accelerate initializer when using exact kernel?

QUESTIONS:
1. The first step: we compute N full line arrangements. That gives an intersection graph. But to get it, we need to intersect 3D planes. Right?
2. How do they handle cases with multiple coplanar polygons?
3. Can we avoid kinetic completely? What if we use labeling instead?
4. What is the exact stop criteria for k intersections? When two polygons intersect at the very beginning, does it count as an intersection? How to count it when we meet an iedge or an ivertex? Can we squeeze through the whole or not?
5. It looks like for open case, we never insert new polygons. Or better, do we need open case at all? Maybe it should be reassigned into closing case! Because, in the original paper, they do not have such a case.
6. Graph-cut based surface extraction does not guarantee a 2-manifold. What if we simply stop when we meet a roof polygon or another wall polygon? In this case, we do not guarantee convex polyhedra but who cares?
7. Do we use the trick from the paper when inserting only the first three events (corresponding to the three closest lines) of an active vertex in the queue instead of all possible events?

ANSWERS:
1. Do we need to intersect all 3D planes in order to find N full line arrangements?
- Yes.

2. How do you handle cases with multiple coplanar polygons?
- Merge polygons, which are at the same cell.

3. Can we avoid kinetic completely?
- Probably yes, but there is a problem of how to find the correct criteria that guarantees the volume convexity. But we can use the trick with tagging faces at least for all initial polygons.

4. Do you have any guarantee on the number of final volumes? Is this number optimal for each k?
- No, it is super difficult.

5. When two polygons intersect at the very beginning, does it count as an intersection? Can we squeeze through the whole or not?
- Both are ok.

TODO:
1. Better future directions, compute them exactly and identify if the lines are parallel using the segment coordinates, unify them for all types of events, the function should take two points with the directions and two fixed points and return the future point and the future direction along the edge, one direction must be the limit direction and one direction must be the cropped direction.
2. Precompute occupied iedges while inserting new events.
3. Better graphcut by inserting information, which faces are originated by the roof and facade points.
4. Adaptive time step, e.g. doubled in case we had three steps and did not meet any event.
5. Precompute as mush as possible stuff, e.g. next events.
6. Fix initialization using the same trick I use for adding missing faces.
7. Put adding missing faces and remove faces in two separate files, they can be reused when working on the subdivision.
8. Add the finalizer class.
9. Make the initializer work with the inexact kernel.
10. Better polygon regularizer using the global approach.
11. Graph cut using normals and graphcut using the LOD.
12. Add timing tests, better tests, accelerate the code as much as possible.
13. Try to merge thin volumes in case we are beyond the tolerance value when traversing the volumes.
14. Add interface to insert custom planes for subdivision instead of uniform subdivision, in this case, we can avoid artifacts in the important parts of the model but still be much faster.
15. Try using non-uniform speed for different polygons.
16. Try to avoid initialization and computing the full intersection graph.
17. Try to avoid randomization.
18. Add unconstrained pvertex to ivertex event.
19. Better region growing maybe using the global optimization.
20. Add 3D global regularization.
21. Try to have as few as possible events.
22. Add free-form reconstruction.
23. Add a way to quickly change the k-intersection criteria.
24. Make the code work both with exact and inexact kernels.
25. Add clustering for input clouds.
26. Add automatic learning input parameters.
27. Make the code work with all edge cases.
28. Add missing walls (exterior and interior), add missing roofs, add missing ground.
29. Create LCC.
30. Improve output such that I could return faces iteratively.
31. Make regularization work with exact kernel.
32. Improve time to compile. Split big files into smaller files.
33. KSR 3 -> data_structure (inc. support planes + intersection graph) -> subdivision -> partitioning -> initializer (inc. polygon_splitter) + propagation (inc. event + event_queue) + finalizer (inc. volume extraction); data_structure -> reconstruction -> (shape detection + shape regularization) + visibility + graphcut + model extraction; data_structure -> k_intersection_stop_condition.
34. Compare the timing of our code with the original code.
35. Merge all collinear vertices along input polygons to avoid handling special cases.
36. Implement the function remove_equal_points() before removing collinear points.
37. When merging coplanar polygons, should we give a priority not to the first one but to the biggest one?