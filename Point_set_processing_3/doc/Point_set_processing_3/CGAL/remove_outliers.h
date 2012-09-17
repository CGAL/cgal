namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::remove_outliers()` deletes a user-specified fraction of outliers from the input point set. More specifically, it sorts the input points in increasing order of average squared distances to the \f$ k\f$ nearest neighbors and computes the points with largest value. 

The outliers detection depends on the `k` parameter, specifically the detection of clusters of outliers. 
The number of neighbors should be higher than the size of clusters of outliers in the point set. 
For datasets with no cluster of outliers, this value can be set to a few rings, e.g. 24. Larger value leads to longer computation times. 
For these reasons, we do not provide any default value for this parameter. 

This method modifies the order of input points so as to pack all remaining points first, and returns and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers. 

Example 
-------------- 

See `remove_outliers_example.cpp`. 

CONVERROR StrangeDoc placement

Removes outliers: computes average squared distance to the K nearest neighbors, and sorts the points in increasing order of average distance. 

This method modifies the order of input points so as to pack all remaining points first, and returns an iterator over the first point to remove (see erase-remove idiom). For this reason it should not be called on sorted containers.
\pre k \f$ >\f$= 2.

Template Parameters
--------------

`InputIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `InputIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

iterator over the first point to remove.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `InputIterator` -\f$ >\f$ `Point_3`. `k`: number of neighbors. `threshold_percent`: percentage of points to remove. `kernel`: geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename Kernel> InputIterator remove_outliers(InputIterator first, InputIterator beyond, PointPMap point_pmap, unsigned int k, double threshold_percent, const Kernel& kernel);

} /* namespace CGAL */

