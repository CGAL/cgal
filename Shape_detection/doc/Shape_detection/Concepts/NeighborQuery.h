/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing`
to access neighbors of an item.

\cgalHasModel
- `CGAL::Shape_detection::Point_set::K_neighbor_query`,
- `CGAL::Shape_detection::Point_set::Sphere_neighbor_query`,
- `CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query`
*/
class NeighborQuery {

public:

  /*!
    fills `neighbors` with the indices of all items, which are connected to the
    item with the index `query_index`.

    `CGAL::Shape_detection::Region_growing` calls this function each time when
    a new query item is selected.
  */
  void operator()(
    const std::size_t query_index,
    std::vector<std::size_t>& neighbors) {

  }
};
