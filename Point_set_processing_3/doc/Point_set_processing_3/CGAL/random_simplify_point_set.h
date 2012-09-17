namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::random_simplify_point_set()` randomly deletes a user-specified fraction of the input points. This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers. 

\sa `CGAL::grid_simplify_point_set`

Example 
-------------- 

See `random_simplification_example.cpp`. 

CONVERROR StrangeDoc placement

Randomly deletes a user-specified fraction of the input points. 

This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers.

Template Parameters
--------------

`ForwardIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `ForwardIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

iterator over the first point to remove.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `ForwardIterator` -\f$ >\f$ `Point_3`. `removed_percentage`: percentage of points to remove. `kernel`: geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> ForwardIterator random_simplify_point_set(ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, double removed_percentage, const Kernel& kernel);

} /* namespace CGAL */

