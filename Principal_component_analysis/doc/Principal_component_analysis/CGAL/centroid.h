namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the centroid of a non-empty set of 2D objects. The tag
is used to specify the dimension to be considered from the
objects. `K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.

The function `centroid` computes the (uniform) center of mass of a set
of 2D or 3D bounded objects. In 2D these objects include points,
segments, triangles, iso rectangles, circles and disks. In 3D these
objects include points, segments, triangles, iso cuboids, spheres,
balls and tetrahedra.

The user can also optionally pass an explicit kernel, in case the
default based on `Kernel_traits` is not sufficient. The default
dimension tag is deduced automatically, although the user can pass a
`tag` specifying the dimension of the objects to be considered for the
centroid computation. For example, the default dimension of a
tetrahedron is 3, but specifying a dimension 0 computes the centroid
of the tetrahedron vertices (3D points), specifying a dimension 1
computes the centroid of the tetrahedron edges (3D segments) and
specifying a dimension 2 computes the centroid of the tetrahedron
facets (3D triangles).

\sa `CGAL::barycenter` 

The value type must be either `K::Point_2`, `K::Segment_2`,
`K::Triangle_2`, `K::Rectangle_2` or `K::Circle_2`. To fit a set of
disks the user must call the function with value type `K::Circle_2`
and with dimension tag of 2. The tag must range between
`CGAL::Dimension_tag<0>` and `CGAL::Dimension_tag<2>`.

\pre first != beyond. 
*/
template < typename InputIterator, typename Tag >
K::Point_2
centroid(InputIterator first, InputIterator beyond, const Tag& t);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the centroid of a non-empty set of 2D objects. The tag
is used to specify the dimension to be considered from the
objects. `K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.

The function `centroid` computes the (uniform) center of mass of a set
of 2D or 3D bounded objects. In 2D these objects include points,
segments, triangles, iso rectangles, circles and disks. In 3D these
objects include points, segments, triangles, iso cuboids, spheres,
balls and tetrahedra.

The user can also optionally pass an explicit kernel, in case the
default based on `Kernel_traits` is not sufficient. The default
dimension tag is deduced automatically, although the user can pass a
`tag` specifying the dimension of the objects to be considered for the
centroid computation. For example, the default dimension of a
tetrahedron is 3, but specifying a dimension 0 computes the centroid
of the tetrahedron vertices (3D points), specifying a dimension 1
computes the centroid of the tetrahedron edges (3D segments) and
specifying a dimension 2 computes the centroid of the tetrahedron
facets (3D triangles).

\sa `CGAL::barycenter` 

The value type must be either `K::Point_2`, `K::Segment_2`,
`K::Triangle_2`, `K::Rectangle_2` or `K::Circle_2`. To fit a set of
disks the user must call the function with value type `K::Circle_2`
and with dimension tag of 2. The tag must range between
`CGAL::Dimension_tag<0>` and `CGAL::Dimension_tag<2>`.

\pre first != beyond. 
*/
template < typename InputIterator, typename K, typename Tag >
K::Point_2
centroid(InputIterator first, InputIterator beyond, const K & k, const Tag& t);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the centroid of a non-empty set of 3D objects. The tag
is used to specify the dimension to be considered from the
objects. `K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.

The function `centroid` computes the (uniform) center of mass of a set
of 2D or 3D bounded objects. In 2D these objects include points,
segments, triangles, iso rectangles, circles and disks. In 3D these
objects include points, segments, triangles, iso cuboids, spheres,
balls and tetrahedra.

The user can also optionally pass an explicit kernel, in case the
default based on `Kernel_traits` is not sufficient. The default
dimension tag is deduced automatically, although the user can pass a
`tag` specifying the dimension of the objects to be considered for the
centroid computation. For example, the default dimension of a
tetrahedron is 3, but specifying a dimension 0 computes the centroid
of the tetrahedron vertices (3D points), specifying a dimension 1
computes the centroid of the tetrahedron edges (3D segments) and
specifying a dimension 2 computes the centroid of the tetrahedron
facets (3D triangles).

\sa `CGAL::barycenter` 

The value type must be either `K::Point_3`, `K::Segment_3`,
`Triangle_3`, `Cuboid_3`, `Sphere_3` or `Tetrahedron_3`. To fit a set
of balls the user must call the function with value type `K::Sphere_3`
and with dimension tag of 3. The tag must range between
`CGAL::Dimension_tag<0>` and `CGAL::Dimension_tag<3>`.

\pre first != beyond. 
*/
template < typename InputIterator, typename Tag >
K::Point_3
centroid(InputIterator first, InputIterator beyond, const Tag& t);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the centroid of a non-empty set of 3D objects. The tag
is used to specify the dimension to be considered from the
objects. `K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.

The function `centroid` computes the (uniform) center of mass of a set of 2D or 3D bounded objects. In 2D these objects include points, segments, triangles, iso rectangles, circles and disks. In 3D these objects include points, segments, triangles, iso cuboids, spheres, balls and tetrahedra. 

The user can also optionally pass an explicit kernel, in case the default based on `Kernel_traits` is not sufficient. The default dimension tag is deduced automatically, although the user can pass a `tag` specifying the dimension of the objects to be considered for the centroid computation. For example, the default dimension of a tetrahedron is 3, but specifying a dimension 0 computes the centroid of the tetrahedron vertices (3D points), specifying a dimension 1 computes the centroid of the tetrahedron edges (3D segments) and specifying a dimension 2 computes the centroid of the tetrahedron facets (3D triangles). 

\sa `CGAL::barycenter` 

The value type must be either `K::Point_3`, `K::Segment_3`,
`Triangle_3`, `Cuboid_3`, `Sphere_3` or `Tetrahedron_3`. To fit a set
of balls the user must call the function with value type `K::Sphere_3`
and with dimension tag of 3. The tag must range between
`CGAL::Dimension_tag<0>` and `CGAL::Dimension_tag<3>`.

\pre first != beyond. 
*/
template < typename InputIterator, typename K, typename Tag >
K::Point_3
centroid(InputIterator first, InputIterator beyond, const K & k, const Tag& t);

} /* namespace CGAL */

