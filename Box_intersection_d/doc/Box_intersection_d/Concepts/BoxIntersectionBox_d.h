/*!
\ingroup PkgBoxIntersectionDConcepts
\cgalConcept

The `BoxIntersectionBox_d` concept is used in the context of the intersection
algorithms for sequences of iso-oriented boxes. These algorithms come
with a default traits class that assumes that the boxes are a model of
this `BoxIntersectionBox_d` concept. This concept defines the access functions to
the dimension, the `id`-number, and the boundaries of the box.

\cgalRefines `Assignable`

\cgalHasModel CGAL::Box_intersection_d::Box_d
\cgalHasModel CGAL::Box_intersection_d::Box_with_handle_d

\sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
\sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>`
\sa `BoxIntersectionTraits_d`

*/

class BoxIntersectionBox_d {
public:

/// \name Types
/// @{

/*!
number type to represent the box
boundaries. Allowed are the built-in types `int`, `unsigned
int`, `float`, and `double`.
*/
typedef unspecified_type NT;

/*!
type for the box `id`-number,
must be a model of the `LessThanComparable` concept.
*/
typedef unspecified_type ID;

/*!
returns the dimension of the box.
*/
static int dimension();

/// @}

/// \name Access Functions
/// @{

/*!
returns the unique `id`-number for the `box`.
*/
ID id() const;

/*!
returns the lower boundary in dimension
`d`, \f$ 0 \leq\f$`d`\f$ < \f$`dimension()`.
*/
NT min_coord( int d) const;

/*!
returns the upper boundary in dimension
`d`, \f$ 0 \leq\f$`d`\f$ < \f$`dimension()`.
*/
NT max_coord( int d) const;

/// @}

}; /* end BoxIntersectionBox_d */
