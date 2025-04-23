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

\todo Add the requirements in the concept `ConformingConstrainedDelaunayTriangulationTraits_3`.

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
  * The line type
  */
  using Line_3 = unspecified_type;

  /*!
  * The point type
  */
  using Point_3 = unspecified_type;

  /*!
  * The vector type
  */
  using Vector_3 = unspecified_type;

  /*!
  * The segment type
  */
  using Segment_3 = unspecified_type;

  /*!
  * The triangle type
  */
  using Triangle_3 = unspecified_type;

  /*!
  * The tetrahedron type
  */
  using Tetrahedron_3 = unspecified_type;

  /*!
   */
  using Construct_cross_product_vector_3 = unspecified_type;

  /*!
   */
  using Construct_vector_3 = unspecified_type;

  /*!
   */
  using Construct_sum_of_vectors_3 = unspecified_type;

  /*!
   */
  using Construct_triangle_3 = unspecified_type;

  /*!
   */
  using Construct_vertex_3 = unspecified_type;

  /*!
   */
  using Orientation_3 = unspecified_type;

  /*!
  * A predicate object that must provide the function operators:

    `bool operator()(Triangle_3 t)`
    `bool operator()(Tetrahedron_3 t)`

    which return true iff the object is degenerate.
   */
  using Is_degenerate_3 = unspecified_type;

  /// @}

  /// \name Operations
  /// The following functions give access to the predicate and construction objects:
  /// @{

  /*!
   */
  Construct_cross_product_vector_3 construct_cross_product_vector_3_object();

  /*!
   */
  Construct_vector_3 construct_vector_3_object();

  /*!
   */
  Construct_sum_of_vectors_3 construct_sum_of_vectors_3_object();

  /*!
   */
  Construct_triangle_3 construct_triangle_3_object();

  /*!
   */
  Construct_vertex_3 construct_vertex_3_object();

  /*!
   */
  Orientation_3 orientation_3_object();

  /*!
   */
  Is_degenerate_3 is_degenerate_3_object();

  /// @}

};
