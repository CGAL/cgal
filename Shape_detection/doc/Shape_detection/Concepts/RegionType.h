/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing`
to maintain a region.

A region is represented by a set of `indices` of the items, which are included in
this region.

\cgalHasModel
- `CGAL::Shape_detection::Point_set::Least_squares_line_fit_region`
- `CGAL::Shape_detection::Point_set::Least_squares_circle_fit_region`
- `CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region`
- `CGAL::Shape_detection::Point_set::Least_squares_sphere_fit_region`
- `CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region`
- `CGAL::Shape_detection::Segment_set::Least_squares_line_fit_region`
- `CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region`
- `CGAL::Shape_detection::Polyline::Least_squares_line_fit_region`
*/
class RegionType {

public:

  /*!
    checks if the item with the index `index_to`, which is a neighbor of the item
    with the index `index_from`, can be added to the region represented by `indices`.

    `CGAL::Shape_detection::Region_growing` calls this function each time when
    trying to add a new item to a region. If this function returns `true`, the
    item with the index `index_to`, is added to the region, otherwise ignored.
  */
  bool is_part_of_region(
    const std::size_t index_from,
    const std::size_t index_to,
    const std::vector<std::size_t>& indices) {

  }

  /*!
    checks if the region represented by `indices` satisfies all necessary conditions.

    `CGAL::Shape_detection::Region_growing` calls this function at the end of each
    propagation phase. If this function returns `true`, the region is accepted,
    otherwise rejected. If the region is rejected, all its items are released and
    available for region growing again.
  */
  bool is_valid_region(
    const std::vector<std::size_t>& indices) {

  }

  /*!
    enables to update any information that is maintained with the region
    represented by `indices`.

    `CGAL::Shape_detection::Region_growing` calls this function each time when a
    new seed item is selected. This case can be identified by checking the
    condition `indices.size() == 1`. This function is also called periodically
    when enlarging the region. This case can be identified by checking the
    condition `indices.size() > 1`.

    This function also returns a Boolean at the first call when a new region
    with one seed item is being created. When it is `true`, the new region is
    further propagated, otherwise, it is rejected.
  */
  bool update(
    const std::vector<std::size_t>& indices) {

  }

  /*!
    enables to update any information that is maintained with the region
    represented by `indices`.

    \deprecated This is deprecated function that was used in all versions before \cgal 5.4.
    The change from `void` in the old version to `bool` in the new version affects
    only user-defined region types.
  */
  void update(
    const std::vector<std::size_t>& indices) {

  }
};
