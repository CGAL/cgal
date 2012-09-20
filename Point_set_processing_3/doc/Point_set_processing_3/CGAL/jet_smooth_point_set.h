namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::jet_smooth_point_set()` smoothes a point set by fitting for each point a jet surface and projecting it onto the jet. The default jet surface is a quadric. 

The output of the smoothing algorithm highly depends on the `k` parameter. 
The number of neighbors controls the size of the point subset considered for jet fitting at each input point. As this parameter is application-specific we do not provide any default value. Larger values lead to smoother point sets and are more time consuming. For point sets with limited noise this value can be set to small number such as 24. For noisy point sets this value must be increased. 

As this method relocates the points, it should not be called on containers sorted w.r.t. point locations.

\pre k \f$ >\f$= 2.

\tparam InputIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `InputIterator`  is convertible to `Point_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\param first, beyond iterator range of the input points. 
\param point_pmap property map `InputIterator` -> `Point_3`. 
\param k number of neighbors. 
\param kernel geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> void jet_smooth_point_set(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, const Kernel& kernel, unsigned int degree_fitting = 2, unsigned int degree_monge = 2);

} /* namespace CGAL */

