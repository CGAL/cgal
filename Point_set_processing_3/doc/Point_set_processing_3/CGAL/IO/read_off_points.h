namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::read_off_points()` reads the point from a .off ASCII stream. More specifically, it reads only the point locations and ignores all point attributes available from the stream. 
`CGAL::read_off_points_and_normals()` reads the points as well as the normals (if available) from a .off ASCII stream. 
In both cases the other primitives (segments, faces) are ignored. 

\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::write_off_points`

Example 
-------------- 

See `read_write_xyz_point_set_example.cpp`. 

CONVERROR StrangeDoc placement

Reads points (position only) from a .off ASCII stream. The function expects for each point a line with the x y z position. If the position is followed by the nx ny nz normal, then the normal will be ignored. Faces are ignored.

Template Parameters
--------------

`OutputIterator`: iterator over output points. `PointPMap`: is a model of `boost::WritablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `OutputIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

true on success.

Parameters
--------------

`stream`: input stream. `output`: output iterator over points. `point_pmap`: property map `OutputIterator` -\f$ >\f$ `Point_3`. `kernel`: geometric traits.

*/
template<typename OutputIterator, typename PointPMap, typename Kernel> bool read_off_points(std::istream& stream, OutputIterator output, PointPMap point_pmap, const Kernel& kernel);

/*!
\ingroup PkgPointSetProcessing

`CGAL::read_off_points()` reads the point from a .off ASCII stream. More specifically, it reads only the point locations and ignores all point attributes available from the stream. 
`CGAL::read_off_points_and_normals()` reads the points as well as the normals (if available) from a .off ASCII stream. 
In both cases the other primitives (segments, faces) are ignored. 

\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::write_off_points`

Example 
-------------- 

See `read_write_xyz_point_set_example.cpp`. 

CONVERROR StrangeDoc placement

Reads points (positions + normals, if available) from a .off ASCII stream. The function expects for each point a line with the x y z position, optionally followed by the nx ny nz normal. Faces are ignored.

Template Parameters
--------------

`OutputIterator`: iterator over output points. `PointPMap`: is a model of `boost::WritablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `OutputIterator` `value_type` is convertible to `Point_3<Kernel>`. `NormalPMap`: is a model of `boost::WritablePropertyMap` with a `value_type` = `Vector_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

true on success.

Parameters
--------------

`stream`: input stream. `output`: output iterator over points. `point_pmap`: property map `OutputIterator` -\f$ >\f$ `Point_3`. `normal_pmap`: property map `OutputIterator` -\f$ >\f$ `Vector_3`. `kernel`: geometric traits.

*/
template<typename OutputIterator, typename PointPMap, typename NormalPMap, typename Kernel> bool read_off_points_and_normals(std::istream& stream, OutputIterator output, PointPMap point_pmap, NormalPMap normal_pmap, const Kernel& kernel);

} /* namespace CGAL */

