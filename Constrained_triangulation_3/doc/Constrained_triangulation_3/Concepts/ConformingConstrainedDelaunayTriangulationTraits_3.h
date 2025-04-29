/*!
\ingroup PkgConstrainedTriangulation3Concepts
\cgalConcept

The concept `ConformingConstrainedDelaunayTriangulationTraits_3` specifies the requirements
for the geometric traits class of the triangulation used as the first template
parameter `Triangulation_3` in the function template
`CGAL::make_conforming_constrained_Delaunay_triangulation_3()`.

\cgalRefines{DelaunayTriangulationTraits_3, ProjectionTraitsGeometricTraits_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Exact_predicates_inexact_constructions_kernel (recommended)}
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

*/
class ConformingConstrainedDelaunayTriangulationTraits_3 {
public:

  /// \name Operations
  /// The following functions give access to the predicate and construction objects:
  /// @{

  /*!
   * returns a function object model of `Kernel::Angle_3`
   */
  unspecified_type angle_3_object();

  /*!
   * returns a function object model of `Kernel::CompareAngle_3`
   */
  unspecified_type compare_angle_3_object();

  /*!
   * returns a function object model of `Kernel::ComputeScalarProduct_3`
   */
  unspecified_type compute_scalar_product_3_object();

  /*!
   * returns a function object model of `Kernel::ComputeSquaredLength_3`
   */
  unspecified_type compute_squared_length_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructCrossProductVector_3`
   */
  unspecified_type construct_cross_product_vector_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructMidpoint_3`
   */
  unspecified_type construct_midpoint_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructVector_3`
   */
  unspecified_type construct_vector_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructScaledVector_3`
   */
  unspecified_type construct_scaled_vector_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructSumOfVectors_3`
   */
  unspecified_type construct_sum_of_vectors_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructTranslatedPoint_3`
   */
  unspecified_type construct_translated_point_3_object();

  /*!
   * returns a function object model of `Kernel::ConstructVertex_3`
   */
  unspecified_type construct_vertex_3_object();

  /*!
   * returns a predicate object model of `Kernel::IsDegenerate_3`
   */
  unspecified_type is_degenerate_3_object();

  /// @}

};
