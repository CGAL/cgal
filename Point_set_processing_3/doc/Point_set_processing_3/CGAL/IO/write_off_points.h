namespace CGAL {

/*!
\ingroup PkgPointSetProcessing

`CGAL::write_off_points()` saves the points of an iterator range to a .off ASCII stream. More specifically, it saves only the point locations and ignores all other attributes. 
`CGAL::write_off_points_and_normals()` saves the points as well the normals of an iterator range to a .off ASCII stream. 

\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`

Example 
-------------- 

See `read_write_xyz_point_set_example.cpp`. 

CONVERROR StrangeDoc placement

Saves the [first, beyond) range of points (positions only) to a .off ASCII stream. The function writes for each point a line with the x y z position.

Template Parameters
--------------

`ForwardIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `ForwardIterator` `value_type` is convertible to `Point_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

true on success.

Parameters
--------------

`stream`: output stream. `first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `ForwardIterator` -\f$ >\f$ `Point_3`. `kernel`: geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename Kernel> bool write_off_points(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, const Kernel& kernel);

/*!
\ingroup PkgPointSetProcessing

`CGAL::write_off_points()` saves the points of an iterator range to a .off ASCII stream. More specifically, it saves only the point locations and ignores all other attributes. 
`CGAL::write_off_points_and_normals()` saves the points as well the normals of an iterator range to a .off ASCII stream. 

\sa `CGAL::read_xyz_points`
\sa `CGAL::write_xyz_points`
\sa `CGAL::read_off_points`

Example 
-------------- 

See `read_write_xyz_point_set_example.cpp`. 

CONVERROR StrangeDoc placement

Saves the [first, beyond) range of points (positions + normals) to a .off ASCII stream. The function writes for each point a line with the x y z position followed by the nx ny nz normal.
\pre normals must be unit vectors.

Template Parameters
--------------

`ForwardIterator`: iterator over input points. `PointPMap`: is a model of `boost::ReadablePropertyMap` with a `value_type` = `Point_3<Kernel>`. It can be omitted if `ForwardIterator` `value_type` is convertible to `Point_3<Kernel>`. `NormalPMap`: is a model of `boost::WritablePropertyMap` with a `value_type` = `Vector_3<Kernel>`. `Kernel`: Geometric traits class. It can be omitted and deduced automatically from `PointPMap` `value_type`.

Returns
--------------

true on success.

Parameters
--------------

`stream`: output stream. `first`: iterator over the first input point. `beyond`: past-the-end iterator over the input points. `point_pmap`: property map `ForwardIterator` -\f$ >\f$ `Point_3`. `normal_pmap`: property map `ForwardIterator` -\f$ >\f$ `Vector_3`. `kernel`: geometric traits.

*/
template<typename ForwardIterator, typename PointPMap, typename NormalPMap, typename Kernel> bool write_off_points_and_normals(std::ostream& stream, ForwardIterator first, ForwardIterator beyond, PointPMap point_pmap, NormalPMap normal_pmap, const Kernel& kernel);

} /* namespace CGAL */

