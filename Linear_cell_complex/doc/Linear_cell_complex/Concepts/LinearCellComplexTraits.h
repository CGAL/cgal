
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

Required types and functors for the `LinearCellComplexTraits` concept. This geometric traits concept is used in the \link CGAL::Linear_cell_complex_for_combinatorial_map `Linear_cell_complex_for_combinatorial_map`\endlink and \link CGAL::Linear_cell_complex_for_generalized_map `Linear_cell_complex_for_generalized_map`\endlink classes.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Linear_cell_complex_traits `CGAL::Linear_cell_complex_traits<d\,K>`\endlink}
\cgalHasModelsEnd

\sa `CGAL::Linear_cell_complex_for_combinatorial_map<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`

*/

class LinearCellComplexTraits {
public:

/// \name Constants
/// @{

/*! The ambient dimension (must be > 1).
*/
static unsigned int ambient_dimension;

/// @}

/// \name Types
/// @{

/*!
a number type that is a model of FieldNumberType.
*/
typedef unspecified_type FT;

/*!
point type.
*/
typedef unspecified_type Point;

/*!
vector type.
*/
typedef unspecified_type Vector;

/// @}

/// \name Constructions
/// @{

/*!
Functor that provides
\link LinearCellComplexTraits::Point `Point`\endlink `operator() (const` \link Point `Point`\endlink`& p, const` \link LinearCellComplexTraits::Vector `Vector`\endlink`& v)`,
which constructs the translation of point `p` by vector `v`, and
\link LinearCellComplexTraits::Point `Point`\endlink `operator() (const CGAL::Origin&, const` \link LinearCellComplexTraits::Vector `Vector`\endlink`& v)`,
which constructs the translation of a point at the origin by vector `v` (used in \link LinearCellComplex::barycenter `barycenter()`\endlink).
*/
typedef unspecified_type Construct_translated_point;

/*!
Functor that provides \link LinearCellComplexTraits::Vector ` Vector `\endlink `operator() (const ` \link Point `Point`\endlink`& p1, const ` \link Point ` Point`\endlink`& p2)`
which constructs a vector as the difference of points `p2-p1`, and
\link LinearCellComplexTraits::Vector ` Vector `\endlink `operator() (const CGAL::Origin&, const ` \link Point ` Point`\endlink`& p)`
which constructs a vector as the difference of point `p` and a point at the origin
(used in \link CGAL::barycenter `barycenter`\endlink
and `CGAL::import_from_plane_graph`).
*/
typedef unspecified_type Construct_vector;

/*!
Functor that provides \link LinearCellComplexTraits::Vector ` Vector `\endlink `operator() (const` \link LinearCellComplexTraits::Vector ` Vector`\endlink`& v1, const` \link LinearCellComplexTraits::Vector ` Vector`\endlink`& v2)`
which constructs a vector as the sum of vectors `v1+v2`
(used in \link CGAL::barycenter `barycenter`\endlink,
`CGAL::compute_normal_of_cell_0`
and `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_sum_of_vectors;

/*!
Functor that provides \link LinearCellComplexTraits::Vector ` Vector `\endlink `operator() (const` \link LinearCellComplexTraits::Vector ` Vector`\endlink`& v, ` \link LinearCellComplexTraits::FT `FT`\endlink `scale)`
which constructs a vector equal to vector `v` scaled by `scale` factor
(used in \link CGAL::barycenter `barycenter`\endlink,
`CGAL::compute_normal_of_cell_0` and `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_scaled_vector;

/*!
Functor that provides \link LinearCellComplexTraits::Point `Point `\endlink `operator() (const ` \link Point `Point`\endlink`& p1, const ` \link Point `Point`\endlink`& p2)`
which constructs the midpoint of points `p1` and `p2`
(used in \link CGAL::barycenter `barycenter`\endlink).
*/
typedef unspecified_type Construct_midpoint;

/// @}

/// \name
/// If `ambient_dimension==2`
/// @{

/*!
a model of \link Kernel::Direction_2 `Direction_2`\endlink.
*/
typedef unspecified_type Direction_2;

/*!
a model of \link Kernel::ConstructDirection_2 `ConstructDirection_2`\endlink (used in `CGAL::import_from_plane_graph`).
*/
typedef unspecified_type Construct_direction_2;

/// @}

/// \name
/// If `ambient_dimension==3`
/// @{

/*!
a model of \link Kernel::ConstructNormal_3 `ConstructNormal_3`\endlink (used in `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_normal_3;

/*!
a model of \link Kernel::Collinear_3 `Collinear_3`\endlink (used in `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Collinear_3;

/// @}

}; /* end LinearCellComplexTraits */
