namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::jet_smooth_point_set()` smooths a point set by fitting for each point a jet surface and projecting it onto the jet. The default jet surface is a quadric. 

The output of the smoothing algorithm highly depends on the `k` parameter. 
The number of neighbors controls the size of the point subset considered for jet fitting at each input point. As this parameter is application-specific we do not provide any default value. Larger values lead to smoother point sets and are more time consuming. For point sets with limited noise this value can be set to small number such as 24. For noisy point sets this value must be increased. 

Example 
-------------- 

See `jet_smoothing_example.cpp`. 

CONVERROR StrangeDoc placement

Smoothes the [first, beyond) range of points using jet fitting on the k nearest neighbors and reprojection onto the jet. As this method relocates the points, it should not be called on containers sorted w.r.t. point locations.
\pre k \f$ >\f$= 2.

Template Parameters
--------------

`InputIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `InputIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `InputIterator` -\f$ >\f$ `Point_3`. `k`: number of neighbors. `kernel`: geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> void jet_smooth_point_set(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, const Kernel& kernel, unsigned int degree_fitting = 2, unsigned int degree_monge = 2);

} /* namespace CGAL */

