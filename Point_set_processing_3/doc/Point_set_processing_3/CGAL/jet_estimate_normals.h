namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::jet_estimate_normals()` estimates normal directions of all points from the input set by fitting jet surfaces over the  `k` nearest neighbors. The default jet surface is a quadric, and the result is an unoriented normal vector for each input point. 

The output of the normal estimation algorithm highly depends on the `k` parameter. 
The number of neighbors controls the size of the point subset considered for jet fitting at each input point. As this parameter is application-specific we do not provide any default value. Larger values lead to smoother normal fields and are more time consuming. For point sets with limited noise this value can be set to small number such as 18. For noisy point sets this value must be increased. 


\tparam InputIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`.
\tparam NormalPMap is a model of `boost::WritablePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\param first, beyond iterator range of the input points. 
\param point_pmap property map `InputIterator` -> `Point_3`. 
\param normal_pmap property map `InputIterator` -> `Vector_3`. 
\param k number of neighbors. 
\param kernel geometric traits.


\sa `CGAL::pca_estimate_normals`
\sa `CGAL::mst_orient_normals`

*/
template<typename InputIterator, typename PointPMap, typename NormalPMap, typename Kernel> void jet_estimate_normals(InputIterator first, InputIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, unsigned int k, const Kernel& kernel, unsigned int degree_fitting = 2);

} /* namespace CGAL */

