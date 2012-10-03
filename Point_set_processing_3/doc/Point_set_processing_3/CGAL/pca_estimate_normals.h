namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::pca_estimate_normals()` estimates normal directions at all points of an input point set by linear least squares fitting of a plane over their `k` nearest neighbors. The result is an unoriented normal for each input point. 

The output of the normal estimation algorithm highly depends on the `k` parameter. 
The number of neighbors controls the size of the point subset considered for plane fitting at each input point. As this parameter is application-specific we do not provide any default value. For noise-free point sets this value can be set to a small number, e.g., 18. Larger values (e.g., 30 or more) lead to smoother normal fields and are more time consuming. We thus recommend using them only for noisy data sets. 

\pre k \f$ >\f$= 2.

\tparam InputIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`. 
\tparam NormalPMap is a model of `boost::WritablePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.


\param first, beyond iterator range of the input points. 
\param point_pmap property map `InputIterator` -> `Point_3`. 
\param normal_pmap property map `InputIterator` -> `Vector_3`. 
\param k number of neighbors. 
\param kernel geometric traits.


\sa `CGAL::jet_estimate_normals`
\sa `CGAL::mst_orient_normals`

*/
template<typename InputIterator, typename PointPMap, typename NormalPMap, typename Kernel> void pca_estimate_normals(InputIterator first, InputIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, unsigned int k, const Kernel& kernel);

} /* namespace CGAL */

