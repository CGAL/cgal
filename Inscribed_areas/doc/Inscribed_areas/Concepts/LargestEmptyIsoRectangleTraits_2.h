
/*!
\ingroup PkgInscribedAreasConcepts
\cgalConcept

The concept `LargestEmptyIsoRectangleTraits_2` describes the set of requirements to be
fulfilled by any class used to instantiate the template parameter of
the class `Largest_empty_iso_rectangle_2<T>`.
This concept provides the types of the geometric primitives used in
this class and some function object types for the required
predicates on those primitives.

\cgalHasModel `CGAL::Cartesian`
\cgalHasModel `CGAL::Homogeneous`

\sa `CGAL::Largest_empty_iso_rectangle_2<Traits>`

*/

class LargestEmptyIsoRectangleTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type.
*/
typedef unspecified_type Point_2;

/*!
The iso rectangle type.
*/
typedef unspecified_type Iso_rectangle_2;

/*!
Predicate object. Must provide
the operator
`Comparison_result operator()(Point_2 p, Point_2 q)`
which returns
`SMALLER, EQUAL` or `LARGER`
according ding to the
\f$ x\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
Predicate object. Must provide
the operator
`Comparison_result operator()(Point_2 p, Point_2 q)`
which returns
`SMALLER, EQUAL` or `LARGER`
according to the
\f$ y\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/*!
Predicate object. Must provide
the operator
`bool operator()(Point_2 p, Point_2 q)`
which returns
whether `p` is less than `q` according to their \f$ x\f$-ordering.
*/
typedef unspecified_type Less_x_2;

/*!
Predicate object. Must provide
the operator
`bool operator()(Point_2 p, Point_2 q)`
which returns
whether `p` is less than `q` according to their \f$ y\f$-ordering.
*/
typedef unspecified_type Less_y_2;

/// @}

/// \name Creation
/// Only a default constructor, copy constructor and an assignement
/// operator are required. Note that further constructors can be
/// provided.
/// @{

/*!
Default constructor.
*/
LargestEmptyIsoRectangleTraits_2();

/*!
Copy constructor
*/
LargestEmptyIsoRectangleTraits_2(LargestEmptyIsoRectangleTraits_2);

/*!
Assignment operator.
*/
LargestEmptyIsoRectangleTraits_2 operator=(LargestEmptyIsoRectangleTraits_2 gtr);

/// @}

/// \name Predicate functions
/// The following functions give access to the predicate and constructor objects.
/// @{

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();

/*!

*/
Less_x_2 less_x_2_object();

/*!

*/
Less_y_2 less_y_2_object();

/// @}

}; /* end LargestEmptyIsoRectangleTraits_2 */

