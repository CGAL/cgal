/*! @brief a signed integral type
 * \relates ConstrainedDelaunayTriangulationVertexBase_3
 */
using CDT_3_face_index = int;

/*!
\cgalConcept

The concept `ConstrainedDelaunayTriangulationVertexBase_3` refines the concept
`TriangulationVertexBase_3` and is the base vertex class for
the `CGAL::Constrained_Delaunay_triangulation_3` class template.

\cgalRefines{TriangulationVertexBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Constrained_Delaunay_triangulation_vertex_base_3}
\cgalHasModelsEnd

\todo add the requirements in the concept `ConstrainedDelaunayTriangulationVertexBase_3`

\sa `ConstrainedDelaunayTriangulationCellBase_3`
*/
class ConstrainedDelaunayTriangulationVertexBase_3 {
public:
  void set_on_constraint(auto constraint_id);

  int number_of_incident_constraints() const;

  void set_mark(CDT_3_vertex_marker marker);

  void clear_marks();

  void clear_mark(CDT_3_vertex_marker marker);

  bool is_marked(CDT_3_vertex_marker marker) const;

  bool is_marked() const;

  template<typename Triangulation>
  auto constraint_id(const Triangulation&) const;

  void set_Steiner_vertex_in_face(CDT_3_face_index face_index);

  CDT_3_face_index face_index() const;

  CDT_3_vertex_type vertex_type() const;
  void set_vertex_type(CDT_3_vertex_type type);
  bool is_Steiner_vertex_on_edge() const;
  bool is_Steiner_vertex_in_face() const;
};

namespace CGAL {
/// \enum CDT_3_vertex_marker
/// @brief Enum type for vertex markers
/// \addtogroup PkgCT_3_enums
enum class CDT_3_vertex_marker {
  CLEAR = 0,     ///< No marker
  REGION_BORDER, ///< On the border of the region
  REGION_INSIDE, ///< Inside the region
  CAVITY,        ///< In the cavity
  CAVITY_ABOVE,  ///< In the cavity above
  CAVITY_BELOW,  ///< In the cavity below
  nb_of_markers  ///< Number of markers
};

/// \enum CDT_3_vertex_type
/// \addtogroup PkgCT_3_enums
enum class CDT_3_vertex_type { FREE, CORNER, STEINER_ON_EDGE, STEINER_IN_FACE };
}
