namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::write_off_points()` saves the points of an iterator range to a .off ASCII stream. More specifically, it saves only the point locations and ignores all other attributes. The function writes for each point a line.

\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator`  is convertible to `Point_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns `true` on success.

\param stream output stream. 
\param first,beyond iterator range over the input points. 
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param kernel geometric traits.


\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`

*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> 
bool 
write_off_points(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, const Kernel& kernel);

/*!
\ingroup PkgPointSetProcessing


`CGAL::write_off_points_and_normals()` saves the points as well the normals of an iterator range to a .off ASCII stream. The function writes for each point a line with the x y z position followed by the nx ny nz normal.


\tparam ForwardIterator iterator over input points. 
\tparam PointPMap is a model of `boost::ReadablePropertyMap` with  value type `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator`  is convertible to `Point_3<Kernel>`. 
\tparam NormalPMap is a model of `boost::WritablePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns `true` on success.

\param stream output stream. 
\param first, beyond iterator range over the input points. 
\param point_pmap property map `ForwardIterator` -> `Point_3`. 
\param normal_pmap property map `ForwardIterator` -> `Vector_3`. 
\param kernel geometric traits.

\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`

*/
template<typename ForwardIterator, typename PointPMap, typename NormalPMap, typename Kernel> 
bool 
write_off_points_and_normals(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, const Kernel& kernel);

} /* namespace CGAL */

