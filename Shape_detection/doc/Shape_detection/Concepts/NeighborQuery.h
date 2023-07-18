/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing`
to access neighbors of an item.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Shape_detection::Point_set::K_neighbor_query}
\cgalHasModels{CGAL::Shape_detection::Point_set::Sphere_neighbor_query}
\cgalHasModels{CGAL::Shape_detection::Polygon_mesh::Polyline_graph}
\cgalHasModels{CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query}
\cgalHasModelsEnd
*/
class NeighborQuery {

public:

  /// The reference type to the elements of the input range, e.g., a const_iterator of the input range.
  typedef unspecified_type Item;

  /*!
    fills `neighbors` with the `Items` of all items, which are connected to the
    `Item` query.

    `CGAL::Shape_detection::Region_growing` calls this function each time when
    a new query item is selected.
  */
  void operator()(
    Item query_index,
    std::vector<Item>& neighbors) {

  }
};
