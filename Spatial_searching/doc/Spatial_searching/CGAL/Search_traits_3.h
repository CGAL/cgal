namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits_3` can be used as a template parameter of the kd tree
and the search classes.

\tparam GeomTraits must be a model of the concept `SearchGeomTraits_3`,
for example `Simple_cartesian<double>` or `Simple_cartesian<Gmpq>`.

\cgalModels{SearchTraits,RangeSearchTraits}

\sa `Search_traits_2<Kernel>`
\sa `Search_traits<NT,Point,CartesianConstIterator,ConstructCartesianConstIterator,Dim>`

*/
template< typename GeomTraits >
class Search_traits_3 {
public:

/// \name Types
/// @{

/*!
Dimension type.
*/
typedef Dimension_tag<3> Dimension;
/*!
Number type.
*/
typedef GeomTraits::FT FT;

/*!
Point type.
*/
typedef GeomTraits::Point_3 Point_d;

/*!
Iso box type.
*/
typedef GeomTraits::Iso_cuboid_3 Iso_box_d;

/*!
Sphere type.
*/
typedef GeomTraits::Sphere_3 Sphere_d;

/*!
An iterator over the %Cartesian
coordinates.
*/
typedef GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_d;

/*!
A functor with
two function operators, which return the begin and past the end iterator for the %Cartesian coordinates.
The functor for begin has as argument a \link Search_traits_3::Point_d `Point_d`\endlink. The functor for the past the end iterator,
has as argument a \link Search_traits_3::Point_d `Point_d`\endlink and an `int`.
*/
typedef GeomTraits::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_d;

/*!
Functor with operator to construct
the iso box from two points.
*/
typedef GeomTraits::Construct_iso_cuboid_3 Construct_iso_box_d;

/*!
Functor with operator to construct
the center of an object of type \link Search_traits_3::Sphere_d `Sphere_d`\endlink.
*/
typedef GeomTraits::Construct_center_3 Construct_center_d;

/*!
Functor with operator to compute
the squared radius of a an object of type \link Search_traits_3::Sphere_d `Sphere_d`\endlink.
*/
typedef GeomTraits::Compute_squared_radius_3 Compute_squared_radius_d;

/*!
Functor with operator to construct
the vertex with lexicographically smallest coordinates of an object of type \link Search_traits_3::Iso_box_d `Iso_box_d`\endlink.
*/
typedef GeomTraits::Construct_min_vertex_3 Construct_min_vertex_d;

/*!
Functor with operator to construct
the vertex with lexicographically largest coordinates of an object of type \link Search_traits_3::Iso_box_d `Iso_box_d`\endlink.
*/
typedef GeomTraits::Construct_max_vertex_3 Construct_max_vertex_d;

/// @}

}; /* end Search_traits_3 */
} /* end namespace CGAL */
