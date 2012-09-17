namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::mst_orient_normals()` orients the normals of a point set using the propagation of a seed orientation through a minimum spanning tree computed over the Riemannian graph \cite cgal:hddms-srup-92. 

The seed is chosen as the top point of the point set. Its normal is oriented towards +Z axis. The success of the orientation algorithm depends on the `k` parameter. 
The number of neighbors controls the number of candidates to propagate the orientation to around each input point. In general the value 18 works well. With smaller values the propagation may be blocked by large gaps in sparse point sets as the graph may be disconnected. Large values cause problems with points scattered over thin objects as the algorithm may incorrectly propagate the orientation from one side of the object to the other. In presence of disconnected clusters of points the algorithm may fail propagating the orientation from one cluster to the others and may only orient the top cluster. 

This method modifies the order of the input points so as to pack all successfully oriented normals first, and returns an iterator over the first point with an unoriented normal (see erase-remove idiom). For this reason it should not be called on sorted containers. 

\sa `CGAL::pca_estimate_normals`
\sa `CGAL::jet_estimate_normals`

Example 
-------------- 

See `pca_estimate_normals_example.cpp`. 

CONVERROR StrangeDoc placement

Orients the normals of the [first, beyond) range of points using the propagation of a seed orientation through a minimum spanning tree of the Riemannian graph \cite cgal:hddms-srup-92. 

This method modifies the order of input points so as to pack all sucessfully oriented points first, and returns an iterator over the first point with an unoriented normal (see erase-remove idiom). For this reason it should not be called on sorted containers.

Preconditions
--------------

Normals must be unit vectors. k \f$ >\f$= 2.

Template Parameters
--------------

`ForwardIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `ForwardIterator` `value_type` is convertible to `Point_3<Kernel>`. `NormalPMap`: is a model of `boost::ReadWritePropertyMap` with a `value_type` = `Vector_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

iterator over the first point with an unoriented normal.

Parameters
--------------

`first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `ForwardIterator` -\f$ >\f$ `Point_3`. `normal_pmap`: property map `ForwardIterator` -\f$ >\f$ `Vector_3`. `k`: number of neighbors. `kernel`: geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename NormalPMap, typename Kernel> ForwardIterator mst_orient_normals(ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, unsigned int k, const Kernel& kernel);

} /* namespace CGAL */

