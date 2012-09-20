namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::remove_outliers()` deletes a user-specified fraction of outliers from the input point set. More specifically, it sorts the input points in increasing order of average squared distances to the `k` nearest neighbors and computes the points with largest value. 

The outliers detection depends on the `k` parameter, specifically the detection of clusters of outliers. 
The number of neighbors should be higher than the size of clusters of outliers in the point set. 
For datasets with no cluster of outliers, this value can be set to a few rings, e.g. 24. Larger value leads to longer computation times. 
For these reasons, we do not provide any default value for this parameter. 

This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers. 



\pre k \f$>\f$= 2.

\tparam InputIterator iterator over input points. 

\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`. 

\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.


\returns iterator over the first point to remove.


\param first, beyond iterator range of the input points. 
\param point_pmap property map `InputIterator` -> `Point_3`. 
\param k number of neighbors. 
\param threshold_percent percentage of points to remove. 
\param kernel geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> InputIterator remove_outliers(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, double threshold_percent, const Kernel& kernel);

} /* namespace CGAL */

