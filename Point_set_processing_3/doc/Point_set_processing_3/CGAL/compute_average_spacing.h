namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::compute_average_spacing()` computes the average spacing of all points from the input set to their `k` nearest neighbors. This value depends on the `k` parameter which can be set to 6 for isotropically sampled surfaces. 


\pre k \f$ >\f$= 2.


\tparam InputIterator iterator over input points. 

\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`. 

\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.


\returns average spacing (scalar).

\param first, beyond iterator range of the input points. 
\param point_pmap property map `InputIterator` -> `Point_3`. 
\param k number of neighbors. 
\param kernel geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> Kernel::FT compute_average_spacing(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, const Kernel& kernel);

} /* namespace CGAL */

