
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`ParameterizerTraits_3` is a concept of parameterization object for a given type
of mesh, `TriangleMesh`, which must be a model of the `FaceGraph` concept.

Creation
--------------

Construction and destruction are undefined.

\cgalHasModel `CGAL::Surface_mesh_parameterization::Parameterizer_traits_3<TriangleMesh>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

*/

class ParameterizerTraits_3
{
public:

/// \name Types
/// @{

/*!

Export the type of mesh to parameterize.

*/
typedef unspecified_type TriangleMesh;

/// @}

/// \name Constants
/// @{

    /// List of errors detected by this package
  enum Error_code {
    OK,                             ///< Success
    ERROR_EMPTY_MESH,               ///< Input mesh is empty
    ERROR_NON_TRIANGULAR_MESH,      ///< Input mesh is not triangular
    ERROR_NO_TOPOLOGICAL_DISC,      ///< Input mesh is not a topological disc
    ERROR_BORDER_TOO_SHORT,         ///< This border parameterization requires a longer border
    ERROR_NON_CONVEX_BORDER,        ///< This parameterization method requires a convex border
    ERROR_CANNOT_SOLVE_LINEAR_SYSTEM,///< Cannot solve linear system
    ERROR_NO_1_TO_1_MAPPING,        ///< Parameterization failed: no one-to-one mapping
    ERROR_WRONG_PARAMETER           ///< A method received an unexpected parameter
  };
/// @}

}; /* end ParameterizerTraits_3 */

