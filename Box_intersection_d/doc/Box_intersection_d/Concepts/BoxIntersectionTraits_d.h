/*!
\ingroup PkgBoxIntersectionDConcepts
\cgalConcept

The `BoxIntersectionTraits_d` concept is used for the intersection algorithms for
sequences of iso-oriented boxes. This concept defines the access
functions to the dimension, the `id`-number, and the boundaries of
the boxes manipulated in these algorithms.

\cgalRefines `Assignable`
\cgalRefines `DefaultConstructible`

\cgalHasModel CGAL::Box_intersection_d::Box_traits_d

\sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink

*/

class BoxIntersectionTraits_d {
public:

/// \name Types
/// @{

/*!
type used for passing box parameters in the
functions below. Since we support in our algorithms
passing the boxes by value as well as passing them as pointers, this
type can be either `const B&`, `B*`, or `const B*`
respectively, where `B` is the actual box type. The difference
to the box handle type lies in the first case where the box handle
would be `B` where this type is `const B&`.
*/
typedef unspecified_type Box_parameter;

/*!
number type to represent the box
boundaries. Allowed are the built-in types `int`, `unsigned
int`, `float`, and `double`.
*/
typedef unspecified_type NT;

/*!
type for the `id`-number,
model of the `LessThanComparable` concept.
*/
typedef unspecified_type ID;

/*!
returns the dimension of the box.
*/
static int dimension();

/*!
returns the unique `id`-number for the `box`.
*/
static ID id(Box_parameter box);

/*!
returns the lower boundary of the `box` in dimension
`d`, \f$ 0 \leq\f$`d`\f$ < \f$`dimension()`.
*/
static NT min_coord( Box_parameter box, int d);

/*!
returns the upper boundary of the `box` in dimension
`d`, \f$ 0 \leq\f$`d`\f$ < \f$`dimension()`.
*/
static NT max_coord( Box_parameter box, int d);

/// @}

}; /* end BoxIntersectionTraits_d */
