
/*!
\ingroup PkgOptimalTransportationReconstruction2Concepts
\cgalConcept

The concept `OptimalTransportationReconstructionTraits_2` describes the requirements 
for the traits class of `CGAL::Optimal_transportation_reconstruction_2`.

\cgalRefines `DelaunayTriangulationTraits_2` 

\cgalHasModel Any model of the `Kernel` concept
\cgalHasModel `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended)

\sa `CGAL::Optimal_transportation_reconstruction_2`

*/

class OptimalTransportationReconstructionTraits_2 {
public:

  /// \name Types 
  /// @{

  /*!
  A coordinate type.
  The type must provide a copy constructor, assignment, comparison
  operators, negation, multiplication, division and allow the
  declaration and initialization with a small integer and double constant
  (cf. requirements for number types).
  An obvious choice would be coordinate type of the point class.
  */
  typedef unspecified_type FT;

  /*!
  The vector type.
  */
  typedef unspecified_type Vector_2;

  /*!
  A function object to construct a `Point_2`.

  Must provides:
  `Point_2 operator()(FT x, FT y)`,
  which constructs a 2D point from its coordinates.
  */
  typedef unspecified_type Construct_point_2;

  /*!
  A function object to construct a `Vector_2`.

  Must provides:
  `Vector_2 operator()(FT x, FT y)`,
  which constructs a 2D vector from its coordinates,
  and `Vector_2 operator()(Point p0, Point p1)`,
  which constructs a 2D vector from 2 points.
  */
  typedef unspecified_type Construct_vector_2;

  /*!
  A function object to construct a `Line_2`.

  Must provides:
  `Line_2 operator()(Point p, Vector v)`,
  which constructs a 2D line from a point and a vector.
  */
  typedef unspecified_type Construct_line_2;

  /*!
  Must provides:
  `Point_2 operator()(Point_2 p, Vector_2 v)`
  that computes the point `p + v`.
  */
  typedef unspecified_type Construct_translated_point_2;

  /*!
  Must provides:
  `Vector_2 operator()(Vector_2 v, FT c)`
  that computes the vector `c * v`.
  */
  typedef unspecified_type Construct_scaled_vector_2;

  /*!
  Must provides:
  `Vector_2 operator()(Vector_2 v1, Vector_2 v2)`
  that computes the vector `v1 + v2`.
  */
  typedef unspecified_type Construct_sum_of_vectors_2;

  /*!
  Must provides:
  `Point_2 operator()(Line_2 l, Point_2 p)`
  that computes the projection of the point `p` on the line `l`.
  */
  typedef unspecified_type Construct_projected_point_2;

  /*!
  Must provides:
  `bool operator()(Line_2 l, Point_2 p)`
  that returns `true` iff `p` lies on `l`.
  */
  typedef unspecified_type Has_on_2;

  /*!
  Must provides:
  `FT operator()(Vector_2 v1, Vector_2 v2)`
  that computes the scalar product between `v1` and `v2`.
  */
  typedef unspecified_type Compute_scalar_product_2;

  /*!
  Must provides:
  `FT operator()(Vector_2 v)`
  that computes the squared length of `v`.
  */
  typedef unspecified_type Compute_squared_length_2;

  /*!
  Must provides:
  `FT operator()(Point p1, Point p2)`
  that computes the squared distance between `p1` and `p2`.
  */
  typedef unspecified_type Compute_squared_distance_2;

  /// @} 

  /// \name Creation 
  /// @{

  /*!
  %Default constructor. 
  */ 
  OptimalTransportationReconstructionTraits_2(); 

  /*!
  Copy constructor. 
  */ 
  OptimalTransportationReconstructionTraits_2 ( 
  const OptimalTransportationReconstructionTraits_2& ); 

  /*!
  Assignment operator.
  */ 
  OptimalTransportationReconstructionTraits_2& operator= 
  (const OptimalTransportationReconstructionTraits_2& ); 

  /// @} 

  /// \name Access to Predicate and Constructors Objects 
  /// @{
  Construct_point_2 construct_point_2_object();
  Construct_vector_2 construct_vector_2_object();
  Construct_vector_2 construct_vector_2_object();
  Construct_line_2 construct_line_2_object();
  Construct_translated_point_2 construct_translated_point_2_object();
  Construct_scaled_vector_2 construct_scaled_vector_2_object();
  Construct_sum_of_vectors_2 construct_sum_of_vectors_2_object();
  Construct_projected_point_2 construct_projected_point_2_object();
  Has_on_2 has_on_2_object();
  Compute_scalar_product_2 compute_scalar_product_2_object();
  Compute_squared_length_2 compute_squared_length_2_object();
  Compute_squared_distance_2 compute_squared_distance_2_object();

  /// @}

}; /* end OptimalTransportationReconstructionTraits_2 */

