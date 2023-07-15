
/*!
\ingroup PkgStreamLines2Concepts
\cgalConcept

The concept `StreamLinesTraits_2` describes the set of requirements for
the template parameter of the class
`CGAL::Regular_grid_2<StreamLinesTraits_2>`.
This concept provides the types handled by the
`CGAL::Stream_lines_2<VectorField_2, Integrator_2>` class.

\cgalHasModelsBegin
\cgalHasModelsBare{The kernels of \cgal are models for this traits class.}
\cgalHasModelsEnd

*/
class StreamLinesTraits_2 {
public:

/// \name Types
/// @{

/*!
The scalar type.
*/
typedef unspecified_type FT;

/*!
The point type.
*/
typedef unspecified_type Point_2;

/*!
The vector type.
*/
typedef unspecified_type Vector_2;

/// @}

}; /* end StreamLinesTraits_2 */

