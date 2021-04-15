The CGAL Open Source Project is pleased to announce the release 5.0 Beta 2
of CGAL, the Computational Geometry Algorithms Library.

CGAL version 5.0 Beta 2 is a public testing release. It should provide
a solid ground to report bugs that need to be tackled before the
release of the final version of CGAL 5.0 in November.

### General changes

- CGAL 5.0 is the first release of CGAL that requires a C++ compiler
  with the support of C++14 or later. The new list of supported
  compilers is:
  - Visual C++ 14.0 (from Visual Studio 2015 Update 3) or later,
  - Gnu g++ 6.3 or later (on Linux or MacOS),
  - LLVM Clang version 8.0 or later (on Linux or MacOS), and
  - Apple Clang compiler versions 7.0.2 and 10.0.1 (on MacOS).
- Since CGAL 4.9, CGAL can be used as a header-only library, with
  dependencies. Since CGAL 5.0, that is now the default, unless
  specified differently in the (optional) CMake configuration.
- The section "Getting Started with CGAL" of the documentation has
  been updated and reorganized.
- The minimal version of Boost is now 1.57.0.


### [Polygonal Surface Reconstruction](https://doc.cgal.org/5.0/Manual/packages.html#PkgPolygonalSurfaceReconstruction) (new package)

 -   This package provides a method for piecewise planar object reconstruction from point clouds. 
     The method takes as input an unordered point set sampled from a piecewise planar object 
     and outputs a compact and watertight surface mesh interpolating the input point set. 
     The method assumes that all necessary major planes are provided (or can be extracted from 
     the input point set using the shape detection method described in Point Set Shape Detection, 
     or any other alternative methods).The method can handle arbitrary piecewise planar objects 
     and is capable of recovering sharp features and is robust to noise and outliers. See also 
     the associated [blog entry](https://www.cgal.org/2019/08/05/Polygonal_surface_reconstruction/).

### [Shape Detection](https://doc.cgal.org/5.0/Manual/packages.html#PkgShapeDetection) (major changes)
 -   **Breaking change:** The concept `ShapeDetectionTraits` has been renamed to [`EfficientRANSACTraits`](https://doc.cgal.org/5.0/Shape_detection/classEfficientRANSACTraits.html).
 -   **Breaking change:** The `Shape_detection_3` namespace has been renamed to [`Shape_detection`](https://doc.cgal.org/5.0/Shape_detection/annotated.html).
 -   Added a new, generic implementation of region growing. This enables for example applying region growing to inputs such as 2D and 3D point sets,
     or models of the [`FaceGraph`](https://doc.cgal.org/5.0/BGL/classFaceGraph.html) concept. Learn more about this new algorithm with this [blog entry](https://www.cgal.org/2019/07/30/Shape_detection/).

### [dD Geometry Kernel](https://doc.cgal.org/5.0/Manual/packages.html#PkgKernelD)
 -   A new exact kernel, [`Epeck_d`](https://doc.cgal.org/5.0/Kernel_d/structCGAL_1_1Epeck__d.html), is now available.

### 2D and 3D Triangulations

-   **Breaking change:** Several deprecated functions and classes have been
    removed. See the full list of breaking changes in the release
    notes.

-   **Breaking change:** The constructor and the `insert()` function of
    `CGAL::Triangulation_2` or `CGAL::Triangulation_3` which take a range
    of points as argument are now guaranteed to insert the points
    following the order of `InputIterator`.  Note that this change only
    affects the base class `CGAL::Triangulation_[23]` and not any
    derived class, such as `CGAL::Delaunay_triangulation_[23]`.


### [Polygon Mesh Processing](https://doc.cgal.org/latest/Manual/packages.html#PkgPolygonMeshProcessing)
 -   Introduced a [wide range of new functions](https://doc.cgal.org/5.0/Polygon_mesh_processing/index.html#title36) 
     related to location of queries on a triangle mesh,
     such as [`CGAL::Polygon_mesh_processing::locate(Point, Mesh)`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__locate__grp.html#gada09bd8740ba69ead9deca597d53cf15).
     The location of a point on a triangle mesh is expressed as the pair of a face and the barycentric
     coordinates of the point in this face, enabling robust manipulation of locations 
     (for example, intersections of two 3D segments living within the same face).
 -   Added the mesh smoothing function [`smooth_mesh()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaa0551d546f6ab2cd9402bea12d8332a3),
     which can be used to improve the quality of triangle elements based on various geometric characteristics.
 -   Added the shape smoothing function [`smooth_shape()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaaa083ec78bcecf351e04d1bbf460b4a2),
     which can be used to smooth the surface of a triangle mesh, using the mean curvature flow to perform noise removal.
     (See also the new entry in the [User Manual](https://doc.cgal.org/5.0/Polygon_mesh_processing/index.html#title8))

### [Point Set Processing](https://doc.cgal.org/latest/Manual/packages.html#PkgPointSetProcessing3)
 -   **Breaking change**: the API using iterators and overloads for optional parameters (deprecated since
     CGAL 4.12) has been removed. The current (and now only) API uses ranges and Named Parameters.

See https://www.cgal.org/2019/10/31/cgal50-beta2/ for a complete list of changes.
