namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::mst_orient_normals()` orients the normals of a point set using the propagation of a seed orientation through a minimum spanning tree computed over the Riemannian graph.

The seed is chosen as the top point of the point set. Its normal is oriented towards +Z axis. The success of the orientation algorithm depends on the `k` parameter. 
The number of neighbors controls the number of candidates to propagate the orientation to around each input point. In general the value 18 works well. With smaller values the propagation may be blocked by large gaps in sparse point sets as the graph may be disconnected. Large values cause problems with points scattered over thin objects as the algorithm may incorrectly propagate the orientation from one side of the object to the other. In presence of disconnected clusters of points the algorithm may fail propagating the orientation from one cluster to the others and may only orient the top cluster. 

This method modifies the order of the input points so as to pack all successfully oriented normals first, and returns an iterator over the first point with an unoriented normal (see erase-remove idiom). For this reason it should not be called on sorted containers. 


\pre Normals must be unit vectors. k \f$ >\f$= 2.

\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator`  is convertible to `Point_3<Kernel>`. 
\tparam NormalPMap is a model of `boost::ReadWritePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns  iterator over the first point with an unoriented normal.

\param first, beyond iterator range of the input points. 
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param normal_pmap property map `ForwardIterator` -> `Vector_3`. 
\param k number of neighbors. 
\param kernel geometric traits.


\sa `CGAL::pca_estimate_normals`
\sa `CGAL::jet_estimate_normals`


*/
template<typename ForwardIterator, typename PointPMap, typename NormalPMap, typename Kernel> ForwardIterator mst_orient_normals(ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, unsigned int k, const Kernel& kernel);

} /* namespace CGAL */

