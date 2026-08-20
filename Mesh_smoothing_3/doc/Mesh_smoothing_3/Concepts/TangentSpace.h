/*!
\ingroup pkgMeshSmoothing3Concepts
\cgalConcept

The concept `TangentSpace` describes a local tangent 
space use for projection on patches or curves.

\sa `CGAL::Mesh_smoothin_3::ConstructTangentSpace`

*/
class TangentSpace {
private:

/// \name Types
/// @{

/*!
defines Point_3 and Vector_3
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
/// The following operation defines the tangent space
/// @{

/*!
returns the origin of the tangent space
*/
Point_3 origin() const;


/*!
returns the vector use to define the space: normal for a plane (Surface) and direction for a line (Curve).
*/
Vector_3 vector() const;

/*!
returns the projection mode on the tangent space
*/
Mesh_smoothing_3::PROJECTION_WEIGHT_MODE projection_mode() const;

/*!
returns weight to use for projection on tangent space. 
Used only if `projection_mode() == PROJECTION_WEIGHT_MODE::CUSTOM`. 
*/
double custom_weight() const;

/// @}



}; /* end TangentSpace */
