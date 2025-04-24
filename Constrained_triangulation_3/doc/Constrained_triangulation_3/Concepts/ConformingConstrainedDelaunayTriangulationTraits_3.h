/*!
\ingroup PkgConstrainedTriangulation3Concepts
\cgalConcept

The concept `ConformingConstrainedDelaunayTriangulationTraits_3` specifies the requirements
for the geometric traits class of the triangulation used as the first template
parameter `Triangulation_3` in the function template
`CGAL::make_conforming_constrained_Delaunay_triangulation_3()`.

\cgalRefines{DelaunayTriangulationTraits_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Exact_predicates_inexact_constructions_kernel (recommended)}
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\todo mention the requirements of `Projection_traits_3<Gt>`

*/
class ConformingConstrainedDelaunayTriangulationTraits_3 {
public:
  /// \name Types
  /// @{

  /*!
  * The number type
  */
  using FT = unspecified_type;

  /*!
  * The vector type
  */
  using Vector_3 = unspecified_type;

  /*!
   * A constructor object model of `ConstructVertex_3`
   */
  using Construct_vertex_3 = unspecified_type;

  /// @}

  /// \name Operations
  /// The following functions give access to the predicate and construction objects:
  /// @{

  /*!
   * returns a function object model of `Angle_3`
   */
  unspecified_type angle_3_object();

  /*!
   * returns a function object model of `CompareAngle_3
   */
  unspecified_type compare_angle_3_object();

  /*!
   * returns a function object model of `ComputeScalarProduct_3`
   */
  unspecified_type compute_scalar_product_3_object();

  /*!
   * returns a function object model of `ComputeSquaredLength_3`
   */
  unspecified_type compute_squared_length_3_object();

  /*!
   * returns a function object model of `ConstructCrossProductVector_3`
   */
  unspecified_type construct_cross_product_vector_3_object();

  /*!
   * returns a function object model of `ConstructMidpoint_3`
   */
  unspecified_type construct_midpoint_3_object();

  /*!
   * returns a function object model of `ConstructVector_3`
   */
  unspecified_type construct_vector_3_object();

  /*!
   * returns a function object model of `ConstructScaledVector_3`
   */
  unspecified_type construct_scaled_vector_3_object();

  /*!
   * returns a function object model of `ConstructSumOfVectors_3`
   */
  unspecified_type construct_sum_of_vectors_3_object();

  /*!
   * returns a function object model of `ConstructTranslatedPoint_3`
   */
  unspecified_type construct_translated_point_3_object();

  /*!
   * returns a function object model of `ConstructVertex_3`
   */
  Construct_vertex_3 construct_vertex_3_object();

  /*!
   * returns a predicate object that must provide the function operators:
    `bool operator()(Triangle_3 t)`
    `bool operator()(Tetrahedron_3 t)`
    which return true iff the object is degenerate.
   */
  unspecified_type is_degenerate_3_object();

  /// @}

};
