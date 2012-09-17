namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::compute_average_spacing()` computes the average spacing of all points from the input set to their \f$ k\f$ nearest neighbors. This value depends on the `k` parameter which can be set to 6 for isotropically sampled surfaces. 

Example 
-------------- 

See `average_spacing_example.cpp`. 

CONVERROR StrangeDoc placement

Computes average spacing from k nearest neighbors.
\pre k \f$ >\f$= 2.

Template Parameters
--------------

`InputIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `InputIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

average spacing (scalar).

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `InputIterator` -\f$ >\f$ `Point_3`. `k`: number of neighbors. `kernel`: geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> Kernel::FT compute_average_spacing(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, const Kernel& kernel);

} /* namespace CGAL */

