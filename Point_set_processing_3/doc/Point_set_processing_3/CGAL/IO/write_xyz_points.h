namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::write_xyz_points()` writes points (positions only) to a .xyz ASCII stream. 


\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns  `true` on success.

\param stream` output stream. 
\param first, beyond iterator range over the input points. 
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param kernel geometric traits.


\sa `CGAL::read_xyz_points`
\sa `CGAL::read_off_points`
\sa `CGAL::write_off_points`

*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> 
bool 
write_xyz_points(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, const Kernel& kernel);

/*!
\ingroup PkgPointSetProcessing

`CGAL::write_xyz_points_and_normals()` saves points (positions + normals) to a .xyz ASCII stream. The function writes for each point a line with the x y z position followed by the nx ny nz normal.

\sa `CGAL::read_xyz_points`
\sa `CGAL::read_off_points`
\sa `CGAL::write_off_points`


\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value typ[e of `ForwardIterator` is convertible to `Point_3<Kernel>`. 
\tparam NormalPMap is a model of `boost::WritablePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns `true` on success.

\pram stream output stream. 
\param first iterator range over the input points.  
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param normal_pmap property map `ForwardIterator` -> `Vector_3`. 
\param kernel geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename NormalPMap, typename Kernel> 
bool 
write_xyz_points_and_normals(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, const Kernel& kernel);

} /* namespace CGAL */

