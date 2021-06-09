* Comment the code.
* Check memory leaking.
* Spell check all files.

* Is it ok license-wise that now this package depends on Mesh_2?
* Why the 2D mesher does not generate boundary points exactly on the polygon boundary but only with the precision 10^-5?
* Should we remove tables? Are not they redundant?

* Merge this package with the natural neighbor coordinates package.
* Improve the solver interface and make it a parameter for the harmonic coordinates class.
* Add a demo with visualization of the coordinate functions.

* Remove Point_2 from the template of free functions and get this type directly from the PointRange typename.
* Put in the classes PointMap first and GeomTraits second.
* Remove the free overload with the traits parameter.