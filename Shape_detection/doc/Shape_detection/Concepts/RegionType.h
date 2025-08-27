/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing`
to maintain a region.

A region is represented by a set items, which are included in this region.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Shape_detection::Point_set::Least_squares_line_fit_region}
\cgalHasModels{CGAL::Shape_detection::Point_set::Least_squares_circle_fit_region}
\cgalHasModels{CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region}
\cgalHasModels{CGAL::Shape_detection::Point_set::Least_squares_sphere_fit_region}
\cgalHasModels{CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region}
\cgalHasModels{CGAL::Shape_detection::Segment_set::Least_squares_line_fit_region}
\cgalHasModels{CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region}
\cgalHasModelsEnd
*/
class RegionType {

public:

  /// The parameters of the primitive covering the region.
  typedef unspecified_type Primitive;

  /// The reference type to the elements of the input range, e.g., a `const_iterator` of the input range.  Must be a model of `Hashable`.
  typedef unspecified_type Item;

  // The Region type is defined by a `vector` of items.
  typedef std::vector<Item> Region;

  /*!
    a model of `ReadWritePropertyMap` whose key type is `Item`
    and value type is `std::size_t`. This map associates item of the input range
    to the index of the region it belongs to.
  */
  typedef unspecified_type Region_index_map;

  /*!
    checks if the `Item` `i` can be added to the `Region` represented by `region`.

    `CGAL::Shape_detection::Region_growing` calls this function each time when
    trying to add a new item to a `Region`. If this function returns `true`, the
    item with the index `i`, is added to the `region`, otherwise ignored.
  */
  bool is_part_of_region(
    const Item i,
    const Region &region) {
  }

  /*!
    checks if `region` satisfies all necessary conditions.

    `CGAL::Shape_detection::Region_growing` calls this function at the end of each
    propagation phase. If this function returns `true`, the `region` is accepted,
    otherwise rejected. If the `region` is rejected, all its items are released and
    available for region growing again.
  */
  bool is_valid_region(
    const Region& region) {
  }

  /*!
    provides the last primitive that has been fitted with the region.
  */

  Primitive primitive() const {
  }

  /*!
    enables to update any information about the region represented by the collection of items `region`.

    `CGAL::Shape_detection::Region_growing` calls this function each time when a
    new seed item is selected. This case can be identified by checking the
    condition `region.size() == 1`. This function is also called periodically
    when enlarging the region. This case can be identified by checking the
    condition `region.size() > 1`.

    This function also returns a boolean at the first call when a new `region`
    with one seed item is being created. When it is `true`, the new `region` is
    further propagated, otherwise, it is rejected.
  */
  bool update(
    const Region& region) {
  }

  }
};
