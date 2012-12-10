
/*!
\ingroup PkgSnapRounding2Concepts
\cgalConcept

The concept `SnapRoundingTraits_2` lists the set of requirements that must be fulfilled by 
an instance of the `Traits` template-parameter of 
the free function \ref CGAL::snap_rounding_2() "CGAL::snap_rounding_2<Traits,InputIterator,OutputContainer>()". 
The list includes the nested types of the geometric primitives used in this class and
some function object types for the required predicates on those primitives.

\cgalRefines `ArrangementTraits_2`

\cgalHasModel `CGAL::Snap_rounding_traits_2<Kernel>` 

\sa `CGAL::Snap_rounding_2<Kernel>` 

\todo check generated documentation page and in particular nested concepts
*/

class SnapRoundingTraits_2 {
public:

/// \name Types 
/// @{

/*! 
The number type. This type must fulfill the requirements on 
`FieldNumberType` 
*/ 
typedef Hidden_type FT; 

/*! 
Models the concept `ArrTraits::Point_2`. 
*/ 
typedef Hidden_type Point_2; 

/*! 
Models the concept `ArrTraits::XMonotoneCurve_2`.
*/ 
typedef Hidden_type Segment_2; 

/*! 
Models the concept `SRTraits_2::IsoRectangle_2`
*/ 
typedef Hidden_type Iso_rectangle_2; 

/// @}

/// \name Functor Types
/// @{

/*! 
Models the concept `SRTraits_2::ConstructVertex_2`.
*/ 
typedef Hidden_type Construct_vertex_2; 

/*! 
Models the concept `SRTraits_2::ConstructSegment_2`.
*/ 
typedef Hidden_type Construct_segment_2; 

/*! 
Models the concept `SRTraits_2::ConstructIsoRectangle_2`. 
*/ 
typedef Hidden_type Construct_iso_rectangle_2; 

/*! 
Models the concept `RealEmbeddableTraits::ToDouble`. The precision of this operation is
of not high significance, as it is only used in the implementation of the
heuristic technique to exploit a cluster of kd-trees rather than just one.
*/ 
typedef Hidden_type To_double; 

/*! 
Models the concept `SRTraits_2::CompareX_2`.
*/ 
typedef Hidden_type Compare_x_2; 

/*! 
Models the concept `SRTraits_2::CompareY_2`. 
*/ 
typedef Hidden_type Compare_y_2; 

/*! 
Models the concept `SRTraits_2::Snap_2`.
*/ 
typedef Hidden_type Snap_2; 

/*! 
Models the concept `SRTraits_2::IntegerGridPoint_2`.
*/ 
typedef Hidden_type Integer_grid_point_2; 

/*! 
Models the concept `SRTraits:MinkowskiSumWithPixel_2`.
*/ 
typedef Hidden_type Minkowski_sum_with_pixel_2; 

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
Construct_iso_rectangle_2 construct_iso_rectangle_2_object(); 

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

/*! 

*/ 
Minkowski_sum_with_pixel_2 minkowski_sum_with_pixel_2_object(); 




/*! 
  Represents an iso rectangle
  \cgalRefines `DefaultConstructible`
  \cgalRefines `CopyConstructible`
  \cgalRefines `Assignable`
  \cgalHasModel `SRTraits_2::Iso_rectangle_2`
*/
class IsoRectangle_2
{};
  
/*!
\cgalRefines `AdaptableBinaryFunction`
\cgalHasModel `SRTraits_2::Construct_vertex_2`
*/
class ConstructVertex_2
{
  public:
  /*!
  returns the source or target of `seg`. If `i` modulo 2 is 0,
  the source is returned, otherwise the target is returned.}  
  */
  Point_2 operator()(Segment_2 seg, int i);
};


/*!
  \cgalRefines `AdaptableBinaryFunction`
  \cgalHasModel `SRTraits_2::Construct_segment_2
*/
class ConstructSegment_2
{
  public:
  /*!
  introduces a segment with source `p` and target `q`. The segment
  is directed from the source towards the target.
  */
  Segment_2 operator()(Point_2 p, Point_2 q);
};


/*!
  \cgalRefines `AdaptableQuaternaryFunction`
  \cgalHasModel `SRTraits_2::Construct_iso_rectangle_2`
*/
class ConstructIsoRectangle_2
{
  public:

  /*!
  introduces an iso-oriented rectangle fo whose minimal `x` coordinate
  is the one of `left`, the maximal `x` coordinate is the one of
  `right`, the minimal `y` coordinate is the one of `bottom`,
  the maximal `y` coordinate is the one of `top`.}
  */
  Iso_rectangle_2 operator()(Point_2 left, Point_2 right,
                             Point_2 bottom, Point_2 top);
};


/*!
  \cgalRefines `AdaptableBinaryFunction`
  \cgalHasModel `SRTraits_2::Compare_x_2`
*/
class CompareX_2
{
  public:
/*!
  returns `SMALLER`, `EQUAL`, or `LARGER` according to the `x`-ordering
  of the points `p` and `q`.
*/
  Comparison_result operator()(Point_2 p, Point_2 q)
};

/*!
  \cgalRefines `AdaptableBinaryFunction`
  \cgalHasModel `SRTraits_2::Compare_y_2`
*/
class CompareY_2
{
  public:
  /*!
    returns `SMALLER`, `EQUAL`, or `LARGER` according to the `y`-ordering
    of the points `p` and `q`.
  */
  Comparison_result operator()(Point_2 p, Point_2 q)
};


/*!
  \cgalRefines `AdaptableQuaternaryFunction`
  \cgalHasModel `SRTraits_2::Snap_2`
*/
class Snap_2
{
  public:
  /*!
    rounds a point to a center of a pixel (unit square) in the grid used by
    the Snap Rounding algorithm. Note that no conversion to an integer grid
    is done yet. `p` is the input point, `pixel_size` is the size of the
    pixel of the grid, and `x` and `y` are the `x` and `y`-coordinates of
    the rounded point respectively.
  */
  void operator()(Point_2 p, FT pixel_size, FT &x, FT &y)
};

/*!
  \cgalRefines `AdaptableBinaryFunction`
  \cgalHasModel `SRTraits_2::Integer_grid_point_2`

*/
class IntegerGridPoint_2
{
  public:
  /*!
    converts coordinates into an integer representation where one unit is equal
    to pixel size. For instance, if a point has the coordinates \f$ (3.7,5.3) \f$
    and the pixel size is \f$ 0.5 \f$, then the new point will have the coordinates
    of \f$ (7,10) \f$. Note, however, that the number type remains the same here,
    although integers are represented. `p` is the converted point and
    `pixel_size` is the size of the pixel of the grid.
  */
  Point_2 operator()(Point_2 p, NT pixel_size)
};

/*!
  \cgalRefines `AdaptableTernaryFunction`
  \cgalHasModel `SRTraits_2::Minkowski_sum_with_pixel_2`
*/
class MinkowskiSumWithPixel_2
{
  public:
  /*!
    returns the vertices of a polygon, which is the Minkowski sum of a segment
    and a square centered at the origin with edge size `pixel edge`.
    `vertices_list` is the list of the vertices of the Minkowski sum
    polygon, `s` is the input segment and `unit_square` is the edge
    size of the pixel.
  */
  void operator()(std::list<Point_2>& vertices_list, Segment_2 s, NT unit_square);
};

/// @}



}; /* end SnapRoundingTraits_2 */

