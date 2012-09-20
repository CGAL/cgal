namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::read_xyz_points()` reads the point locations from a .xyz ASCII stream. The function expects for each point a line with the x y z position. If the position is followed by the nx ny nz normal, then the normal will be ignored. The first line may contain the number of points in the file. Empty lines and comments starting by # character are allowed.

`CGAL::read_xyz_points_and_normals()` reads the points as well as the normals (if available) from a .xyz ASCII stream. 

\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`
\sa `CGAL::write_off_points`

\tparam OutputIterator iterator over output points. 
\tparam PointPMap is a model of `boost::WritablePropertyMap` with value type `Point_3<Kernel>`. It can be omitted if the value type of `OutputIterator`  is convertible to `Point_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns `true` on success.

\param stream input stream. 
\param output output iterator over points. 
\param point_pmap property map `OutputIterator` -> `Point_3`. 
\param kernel geometric traits.

*/
template<typename OutputIterator, typename PointPMap, typename Kernel> 
bool 
read_xyz_points(std::istream& stream, OutputIterator output, PointPMap point_pmap, const Kernel& kernel);

/*!
\ingroup PkgPointSetProcessing

`CGAL::read_xyz_points_and_normals()` reads the points as well as the normals (if available) from a .xyz ASCII stream. 
The function expects for each point a line with the x y z position, optionally followed by the nx ny nz normal. The first line may contain the number of points in the file. Empty lines and comments starting by # character are allowed.

\tparam OutputIterator iterator over output points. 
\tparam PointPMap is a model of `boost::WritablePropertyMap` with  value type `Point_3<Kernel>`. It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`. 
\tparam NormalPMap is a model of `boost::WritablePropertyMap` with value type `Vector_3<Kernel>`. 
\tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the value type of `PointPMap`.

\returns `true` on success.

\param stream input stream. 
\param output output iterator over points. 
\param point_pmap property map `OutputIterator` -> `Point_3`. 
\param normal_pmap property map `OutputIterator` -> `Vector_3`. 
\param kernel geometric traits.


\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`
\sa `CGAL::write_off_points`

*/
template<typename OutputIterator, typename PointPMap, typename NormalPMap, typename Kernel> 
bool 
read_xyz_points_and_normals(std::istream& stream, OutputIterator output, PointPMap point_pmap, NormalPMap normal_pmap, const Kernel& kernel);

} /* namespace CGAL */

