namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::random_simplify_point_set()` randomly deletes a user-specified fraction of the input points. This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers. 


\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns iterator over the first point to remove.

\param first, beyond iterator range of the input points. 
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param removed_percentage percentage of points to remove. 
\param kernel geometric traits.

\sa `CGAL::grid_simplify_point_set`


*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> ForwardIterator random_simplify_point_set(ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, double removed_percentage, const Kernel& kernel);

} /* namespace CGAL */

