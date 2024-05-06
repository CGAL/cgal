/*!
\cgalConcept

The concept `ConstrainedDelaunayTriangulationCellBase_3` refines the concept
`TriangulationCellBase_3` and is the base cell class for
the `CGAL::Constrained_Delaunay_triangulation_3` class template.

\cgalRefines{TriangulationCellBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Constrained_Delaunay_triangulation_cell_base_3}
\cgalHasModelsEnd

\todo add the requirements in the concept `ConstrainedDelaunayTriangulationCellBase_3`

\sa `ConstrainedDelaunayTriangulationVertexBase_3`
*/
class ConstrainedDelaunayTriangulationCellBase_3 {
public:
  bool is_marked() const;
  bool is_marked(CGAL::CDT_3_cell_marker m);
  void set_mark(CGAL::CDT_3_cell_marker m);
  void clear_mark(CGAL::CDT_3_cell_marker m);
  void clear_marks();

  bool is_facet_constrained(int i);

  template <typename Facet_handle>
  void set_facet_constraint(int i, CGAL::CDT_3_face_index face_id,
                            Facet_handle facet_2d);

  CGAL::CDT_3_face_index face_constraint_index(int i);

  template <typename CDT_2>
  auto face_2 (const CDT_2& cdt, int i);
};

namespace CGAL {
/// @brief Enum type for cell markers.
/// \ingroup PkgCT_3_enums
enum class CDT_3_cell_marker {
  CLEAR = 0,
  IN_REGION = 1,
  VISITED = 1,
  ON_REGION_BOUNDARY = 2,
  nb_of_markers
};

}