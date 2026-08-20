/*!
\ingroup pkgMeshSmoothing3Concepts
\cgalConcept

The concept `TangentSpace` describes a local tangent 
space use for projection on patches or curves.

\sa `ConstructTangentSpace`

*/
class TangentSpace {
public:

/// \name Types
/// @{

/*!
Geometric traits defining `Point_3` and `Vector_3`.
*/
using Geom_traits = unspecified_type;

/*!
Point type.
*/
using Point_3 = Geom_traits::Point_3;


/*!
Vector type.
*/
using Vector_3 = Geom_traits::Vector_3;

/// @}


/// \name Operations
/// @{

/*!
returns the origin of the tangent space
*/
Point_3 origin() const;


/*!
Returns the vector defining the tangent space:
a normal for a surface and a tangent direction for a curve.
*/
Vector_3 vector() const;

/*!
returns the projection weighting mode.
*/
CGAL::Mesh_smoothing_3::Projection_weight_mode  projection_mode() const;

/*!
returns the custom projection weight.

This value is used only when `projection_mode()` returns
`CGAL::Mesh_smoothing_3::Projection_weight_mode::CUSTOM`.
*/
double custom_weight() const;

/// @}



}; /* end TangentSpace */
