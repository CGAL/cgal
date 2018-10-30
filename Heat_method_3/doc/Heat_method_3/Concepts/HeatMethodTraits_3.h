/*!
\ingroup PkgHeatMethodConcepts

\cgalConcept

The concept `HeatMethodTraits_3` describes the types,
predicates, and constructions required by the traits class parameter of
`CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3`.

\cgalHasModel All the \cgal kernels


*/

class HeatMethodTraits_3
{
public:

/// \name Types
/// @{

  /// The 3D point type.
  typedef unspecified_type Point_3;
  
  /// The 3D vector type.
  typedef unspecified_type Vector_3;
  
  /// The 2D point type.
  typedef unspecified_type Point_2;

  typedef unspecified_type Compute_x_3;
  typedef unspecified_type Compute_y_3;
  typedef unspecified_type Compute_z_3;
  typedef unspecified_type Construct_vector_3;
    typedef unspecified_type Construct_sum_of_vectors_3;
  typedef unspecified_type Construct_scaled_vector_3;
  typedef unspecified_type Construct_cross_product_vector_3;
  typedef unspecified_type Compute_scalar_product_3;
  typedef unspecified_type Compute_squared_distance_3;
  
/// @}


};
