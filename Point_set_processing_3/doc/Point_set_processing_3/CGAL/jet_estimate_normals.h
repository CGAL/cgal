namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::jet_estimate_normals()` estimates normal directions of all points from the input set by fitting jet surfaces over the \f$ k\f$ nearest neighbors. The default jet surface is a quadric, and the result is an unoriented normal vector for each input point. 

The output of the normal estimation algorithm highly depends on the `k` parameter. 
The number of neighbors controls the size of the point subset considered for jet fitting at each input point. As this parameter is application-specific we do not provide any default value. Larger values lead to smoother normal fields and are more time consuming. For point sets with limited noise this value can be set to small number such as 18. For noisy point sets this value must be increased. 

\sa `CGAL::pca_estimate_normals`
\sa `CGAL::mst_orient_normals`

Example 
-------------- 

See `normal_estimation.cpp` example. 

CONVERROR StrangeDoc placement

Estimates normal directions of the [first, beyond) range of points using jet fitting on the k nearest neighbors. The output normals are randomly oriented.
\pre k \f$ >\f$= 2.

Template Parameters
--------------

`InputIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `InputIterator` `value_type` is convertible to `Point_3<Kernel>`. `NormalPMap`: is a model of `boost::WritablePropertyMap` with a `value_type` = `Vector_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `InputIterator` -\f$ >\f$ `Point_3`. `normal_pmap`: property map `InputIterator` -\f$ >\f$ `Vector_3`. `k`: number of neighbors. `kernel`: geometric traits.

*/
template<typename InputIterator, typename PointPMap, typename NormalPMap, typename Kernel> void jet_estimate_normals(InputIterator first, InputIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, unsigned int k, const Kernel& kernel, unsigned int degree_fitting = 2);

} /* namespace CGAL */

