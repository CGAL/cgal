/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 *
 * The concept `HotPixelSnapRoundingTraits_2` lists the set of requirements that must be fulfilled by
 * an instance of the `geom_traits` parameter of
 * the free function `CGAL::hot_pixel_snap_rounding_2()`.
 * The list includes the nested types of the geometric primitives used in this class and
 * some function object types for the required predicates on those primitives.
 *
 * \cgalRefines{AosTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>}
 * \cgalHasModelsEnd
 */

class HotPixelSnapRoundingTraits_2 {
public:

/// \name Types
/// @{

/*! The number type. This type must fulfill the requirements on
 * `FieldNumberType` and `RealEmbeddable`.
 */
typedef unspecified_type FT;

/*! Models the concept `AosTraits_2::Point_2`.
 */
typedef unspecified_type Point_2;

/*! Models the concept `AosTraits_2::XMonotoneCurve_2`.
 */
typedef unspecified_type Segment_2;

/// @}

/// \name Functor Types
/// @{

/*! Models the concept `HPSRTraits_2::ConstructVertex_2`.
 */
typedef unspecified_type Construct_vertex_2;

/*! Models the concept `HPSRTraits_2::ConstructSegment_2`.
 */
typedef unspecified_type Construct_segment_2;

/*! Models the concept `HPSRTraits_2::CompareX_2`.
 */
typedef unspecified_type Compare_x_2;

/*! Models the concept `HPSRTraits_2::CompareY_2`.
 */
typedef unspecified_type Compare_y_2;

/*! Models the concept `HPSRTraits_2::Snap_2`.
 */
typedef unspecified_type Snap_2;

/*! Models the concept `HPSRTraits_2::IntegerGridPoint_2`.
 */
typedef unspecified_type Integer_grid_point_2;

/*! Models the concept `AosTraits_2::ConstructMinVertex_2`.
 */
typedef unspecified_type Construct_min_vertex_2;

/*! Models the concept `AosTraits_2::ConstructMaxVertex_2`.
 */
typedef unspecified_type Construct_max_vertex_2;

/// @}

/// \name Accessing Functor Objects
/// @{

/*!
 */
Construct_vertex_2 construct_vertex_2_object();

/*!
 */
Construct_segment_2 construct_segment_2_object();

/*!
 */
Compare_x_2 compare_x_2_object();

/*!
 */
Compare_y_2 compare_y_2_object();

/*!
 */
Snap_2 snap_2_object();

/*!
 */
Integer_grid_point_2 integer_grid_point_2_object();


/// @}

}; /* end HotPixelSnapRoundingTraits_2 */

namespace HPSRTraits_2 {

/*!
 * \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableBinaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Construct_vertex_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Construct_vertex_2` \endlink}
 *  \cgalHasModelsEnd
 */
class ConstructVertex_2 {
public:

  /*! returns the source or target of `seg`. If `i` modulo 2 is 0,
   * the source is returned, otherwise the target is returned.
   */
  Point_2 operator()(Segment_2 seg, int i);
};

/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableBinaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Construct_segment_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Construct_segment_2`  \endlink}
 * \cgalHasModelsEnd
 */
class ConstructSegment_2 {
public:

  /*! introduces a segment with source `p` and target `q`. The segment
   * is directed from the source towards the target.
   */
  Segment_2 operator()(Point_2 p, Point_2 q);
};

/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableBinaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Compare_x_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Compare_x_2` \endlink}
 * \cgalHasModelsEnd
 */
class CompareX_2 {
public:

  /*! returns `SMALLER`, `EQUAL`, or `LARGER` according to the `x`-ordering
   * of the points `p` and `q`.
   */
  Comparison_result operator()(Point_2 p, Point_2 q)
};

/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableBinaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Compare_y_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Compare_y_2` \endlink}
 * \cgalHasModelsEnd
 */
class CompareY_2 {
public:

  /*! returns `SMALLER`, `EQUAL`, or `LARGER` according to the `y`-ordering
   * of the points `p` and `q`.
   */
  Comparison_result operator()(Point_2 p, Point_2 q)
};

/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableQuaternaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Snap_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Snap_2` \endlink}
 * \cgalHasModelsEnd
 */
class Snap_2 {
public:

  /*! rounds a point to a center of a pixel (unit square) in the grid used by
   * the Snap Rounding algorithm. Note that no conversion to an integer grid
   * is done yet. `p` is the input point, `pixel_size` is the size of the
   * pixel of the grid, and `x` and `y` are the `x` and `y`-coordinates of
   * the rounded point respectively.
   */
  void operator()(Point_2 p, FT pixel_size, FT& x, FT& y)
};

/*! \ingroup PkgSnapRounding2Concepts
 * \cgalConcept
 * \cgalRefines{AdaptableBinaryFunction}
 * \cgalHasModelsBegin
 * \cgalHasModelsBare{\link HotPixelSnapRoundingTraits_2::Integer_grid_point_2 `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>::Integer_grid_point_2` \endlink}
 * \cgalHasModelsEnd
 */
class IntegerGridPoint_2 {
public:

  /*! converts coordinates into an integer representation where one unit is equal
   * to pixel size. For instance, if a point has the coordinates \f$ (3.7,5.3) \f$
   * and the pixel size is \f$ 0.5 \f$, then the new point will have the coordinates
   * of \f$ (7,10) \f$. Note, however, that the number type remains the same here,
   * although integers are represented. `p` is the converted point and
   * `pixel_size` is the size of the pixel of the grid.
   */
  Point_2 operator()(Point_2 p, NT pixel_size)
};

} // end of namespace HPSRTraits_2
