
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `FixedWeightedAlphaShapeTraits_3` describes the requirements
for the geometric traits class of the underlying regular triangulation of a weighted alpha shape with fixed alpha value.

\cgalRefines{RegularTriangulationTraits_3 if the underlying triangulation of the alpha shape is a regular triangulation,
  Periodic_3RegularTriangulationTraits_3 if the underlying triangulation of the alpha shape is a periodic regular triangulation}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended kernel)
*/

class FixedWeightedAlphaShapeTraits_3 {
public:

/// \name Types
/// @{

/*!
`CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>`
*/
typedef unspecified_type Comparison_result;

/*!
An object constructor able to compare the squared radius of the
smallest sphere orthogonal to either four, three, two or one weighted point(s)
to a given value of alpha.
It provides:
- `Comparison_result operator()(Weighted_point_3 ,Weighted_point_3 ,Weighted_point_3 ,Weighted_point_3 )`
- `Comparison_result operator()(Weighted_point_3 ,Weighted_point_3 ,Weighted_point_3 )`
- `Comparison_result operator()(Weighted_point_3 ,Weighted_point_3 )`
- `Comparison_result operator()(Weighted_point_3 )`

*/
typedef unspecified_type Compare_weighted_squared_radius_3;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
FixedWeightedAlphaShapeTraits_3();

/// @}

/// \name Access Functions
/// @{

/*!

*/
Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object();

/// @}

}; /* end FixedWeightedAlphaShapeTraits_3 */

