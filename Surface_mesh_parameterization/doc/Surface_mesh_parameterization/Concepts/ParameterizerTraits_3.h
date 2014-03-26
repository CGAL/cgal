
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`ParameterizerTraits_3` is a concept of parameterization object for a given type of mesh, `Adaptor`, which is a model of the `ParameterizationMesh_3` concept. 

Design Pattern 
-------------- 

`ParameterizerTraits_3` models are Strategies \cgalCite{cgal:ghjv-dpero-95} : they implement a strategy of surface parameterization for models of `ParameterizationMesh_3`. 

Creation 
-------------- 

Construction and destruction are undefined. 

\cgalHasModel `CGAL::Parameterizer_traits_3<ParameterizationMesh_3>`
\cgalHasModel `CGAL::Fixed_border_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Barycentric_mapping_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Discrete_authalic_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Discrete_conformal_map_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::LSCM_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\cgalHasModel `CGAL::Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

\sa `ParameterizationMesh_3`

*/

class ParameterizerTraits_3 {
public:

/// \name Types 
/// @{

/*!

Export the type of mesh to parameterize. 

*/ 
typedef unspecified_type Adaptor; 

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
    ERROR_OUT_OF_MEMORY,            ///< Not enough memory
    ERROR_WRONG_PARAMETER           ///< A method received an unexpected parameter
    };
/// @} 

/// \name Operations 
/// @{

/*!

Compute a one-to-one mapping from a triangular 3D surface `mesh` to a piece of the 2D space. The mapping is linear by pieces (linear in each triangle). The result is the (u, v) pair image of each vertex of the 3D surface. 

\pre `mesh` must be a surface with one connected component and no hole. `mesh` must be a triangular mesh. 

*/ 
Error_code parameterize(Adaptor& mesh); 

/// @}

}; /* end ParameterizerTraits_3 */

