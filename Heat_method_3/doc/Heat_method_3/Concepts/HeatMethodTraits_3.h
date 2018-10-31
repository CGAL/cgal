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

  // The number type. It must be `CopyConstructible` and `DefaultConstructible`,
  // and be constructible from `double`.
  typdef unspecified_type FT;
  
  /// The 3D point type. It must be `CopyConstructible` and `DefaultConstructible`,
  /// and have a constructor with three parameters of a type constructibe from `double`.
  typedef unspecified_type Point_3;
  
  /// The 3D vector type.  It must be `CopyConstructible` and `DefaultConstructible`,
  /// and have a constructor with three parameters of a type constructibe from `double`.
  typedef unspecified_type Vector_3;

  /// Functor with operator: `Vector_3 operator()(Point_3 from, Point_3 to) const`.
  typedef unspecified_type Construct_vector_3;

  ///  Functor with operator: `Vector_3 operator()(Vector_3 v, Vector_3 w) const`.
  typedef unspecified_type Construct_sum_of_vectors_3;
  
  ///  Functor with operator: `Vector_3 operator()(Vector_3 v, double d) const`.
  typedef unspecified_type Construct_scaled_vector_3;
  
  ///  Functor with operator: `Vector_3 operator()(Vector_3 v, Vector_3 w) const`.
  typedef unspecified_type Construct_cross_product_vector_3;
  
  ///  Functor with operator: `double operator()(Vector_3 v, Vector_3 w) const`.
  typedef unspecified_type Compute_scalar_product_3;
  
  ///  Functor with operator: `double operator()(Point_3 p, Point_3 q) const`.
  typedef unspecified_type Compute_squared_distance_3;
  
  
/// @}

/*! \name Operations 
For each of the above function object types, 
`Func_obj_type`, a function must exist with the name 
`func_obj_type_object` that creates an instance of the function or 
predicate object type. For example: 
*/
/// @{

/*!

*/ 
Construct_vector_3 construct_vector_3_object(); 

/// @}

};
