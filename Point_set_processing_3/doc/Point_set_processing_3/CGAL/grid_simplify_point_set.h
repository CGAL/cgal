namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

Function `CGAL::grid_simplify_point_set()` considers a regular grid covering the bounding box of the input point set, and clusters all points sharing the same cell of the grid by picking as representant one arbitrarily chosen point. 

This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers. 

\sa `CGAL::random_simplify_point_set`

Example 
-------------- 

See `grid_simplification_example.cpp`. 

CONVERROR StrangeDoc placement

Merges points which belong to the same cell of a grid of cell size = epsilon. 

This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers.
\pre epsilon \f$ >\f$ 0.

Template Parameters
--------------

`ForwardIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `ForwardIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

iterator over the first point to remove.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `ForwardIterator` -\f$ >\f$ `Point_3`. `epsilon`: tolerance value when merging 3D points. `kernel`: geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> ForwardIterator grid_simplify_point_set(ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, double epsilon, const Kernel& kernel);

} /* namespace CGAL */

