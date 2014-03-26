namespace CGAL {

/*!
\addtogroup PkgPrincipalComponentAnalysisDCentroid

The function `centroid()` computes the (uniform) center of mass of a set
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


\sa \link PkgPrincipalComponentAnalysisDBary `CGAL::barycenter()` \endlink
\sa \link centroid_grp `CGAL::centroid() (Linear Kernel)` \endlink

*/
/// @{

/*!
computes the centroid of a non-empty set of 2D or 3D objects. The tag
is used to specify the dimension to be considered from the
objects. 

\pre first != beyond. 

\returns The return type is either `K::Point_2` or `K::Point_3`,
depending on the dimension of the input objects, where `K` is
\code
CGAL::Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel
\endcode

\cgalHeading{Two Dimensional Input}

The value type must be either `K::Point_2`, `K::Segment_2`,
`K::Triangle_2`, `K::Rectangle_2` or `K::Circle_2`. To fit a set of
disks the user must call the function with value type `K::Circle_2`
and with dimension tag of 2. The tag must range between
`Dimension_tag<0>` and `Dimension_tag<2>`.

\cgalHeading{Three Dimensional Input}

The value type must be either `K::Point_3`, `K::Segment_3`,
`K::Triangle_3`, `K::Cuboid_3`, `K::Sphere_3` or `K::Tetrahedron_3`. To fit a set
of balls the user must call the function with value type `K::Sphere_3`
and with dimension tag of 3. The tag must range between
`Dimension_tag<0>` and `Dimension_tag<3>`.

*/
template < typename InputIterator, typename Tag >
Deduced
centroid(InputIterator first, InputIterator beyond, const Tag& t);

/*!  computes the centroid of a non-empty set of 2D or 3D objects. The
tag is used to specify the dimension to be considered from the
objects.

\pre first != beyond. 

\returns The return type is either `K::Point_2` or `K::Point_3`,
depending on the dimension of the input objects.

\cgalHeading{Two Dimensional Input}

The value type must be either `K::Point_2`, `K::Segment_2`,
`K::Triangle_2`, `K::Rectangle_2` or `K::Circle_2`. To fit a set of
disks the user must call the function with value type `K::Circle_2`
and with dimension tag of 2. The tag must range between
`Dimension_tag<0>` and `Dimension_tag<2>`.

\cgalHeading{Three Dimensional Input}

The value type must be either `K::Point_3`, `K::Segment_3`,
`K::Triangle_3`, `K::Cuboid_3`, `K::Sphere_3` or `K::Tetrahedron_3`. To fit a set
of balls the user must call the function with value type `K::Sphere_3`
and with dimension tag of 3. The tag must range between
`Dimension_tag<0>` and `Dimension_tag<3>`.

*/
template < typename InputIterator, typename K, typename Tag >
Deduced
centroid(InputIterator first, InputIterator beyond, const K & k, const Tag& t);

/// @}

} /* namespace CGAL */

