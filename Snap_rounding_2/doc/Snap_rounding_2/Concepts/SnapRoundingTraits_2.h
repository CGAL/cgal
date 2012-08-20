
/*!
\ingroup PkgSnapRounding2Concepts
\cgalconcept

The concept `SnapRoundingTraits_2` lists the set of requirements that must be fulfilled by 
an instance of the `Traits` template-parameter of 
the function \ref CGAL::snap_rounding_2 "CGAL::snap_rounding_2<Traits,InputIterator,OutputContainer>()". 
This concept provides the types of the geometric primitives used in 
this class and some function object types for the required 
predicates on those primitives. 

\refines ::DefaultConstructible
\refines ::Assignable
\refines ::CopyConstructible
\refines ::SweepLineTraits_2
\refines An instance of this concept is used as the traits class for the `CGAL::Sweep_line_2::get_intersection_points()` operation. The requirements listed below are induced by components of the `CGAL::snap_rounding_2` function other than the call to `CGAL::Sweep_line_2::get_intersection_points()`. Naturally, some of them may already be listed in `SweepLineTraits_2`. 

\hasModel CGAL::Snap_rounding_traits 

\sa `CGAL::Snap_rounding_2<Traits>` 

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
The point type. 
*/ 
typedef Hidden_type Point_2; 

/*! 
The segment type. 
*/ 
typedef Hidden_type Segment_2; 

/*! 
The iso-rectangle type. 
*/ 
typedef Hidden_type Iso_rectangle_2; 

/*! 
Function object. Must provide the operator 
`Point_2 operator()(Segment_2 seg, int i)`, which returns the source or 
target of `seg`. If `i` modulo 2 is 0, the source is returned, 
otherwise the target is returned. 
*/ 
typedef Hidden_type Construct_vertex_2; 

/*! 
Function object. Must provide the operator 
`Segment_2 operator()(Point_2 p, Point_2 q)`, which introduces a segment 
with source `p` and target `q`. The segment is directed from the 
source towards the target. 
*/ 
typedef Hidden_type Construct_segment_2; 

/*! 
Function object. Must provide the 
operator 
`Iso_rectangle_2 operator()(Point_2 left, Point_2 right, Point_2 bottom, 
Point_2 top)`, which introduces an iso-oriented rectangle fo whose minimal 
\f$ x\f$ coordinate is the one of `left`, the maximal \f$ x\f$ coordinate is the one 
of `right`, the minimal \f$ y\f$ coordinate is the one of `bottom`, the 
maximal \f$ y\f$ coordinate is the one of `top`. 
*/ 
typedef Hidden_type Construct_iso_rectangle_2; 

/*! 
Function object. Must provide the operator 
`double operator()(FT)`, which computes an approximation of a given number 
of type `FT`. The precision of this operation is of not high significance, 
as it is only used in the implementation of the heuristic technique to exploit 
a cluster of kd-trees rather than just one. 
*/ 
typedef Hidden_type To_double; 

/*! 
Function object. Must provide the operator 
`Comparison_result operator()(Point_2 p, Point_2 q)` 
which returns 
`SMALLER, EQUAL` or `LARGER` according to the 
\f$ x\f$-ordering of points `p` and `q`. 
*/ 
typedef Hidden_type Compare_x_2; 

/*! 
Function object. Must provide the operator 
`Comparison_result operator()(Point_2 p, Point_2 q)` 
which returns 
`SMALLER, EQUAL` or ` LARGER` 
according to the 
\f$ y\f$-ordering of points `p` and `q`. 
*/ 
typedef Hidden_type Compare_y_2; 

/*! 
Rounds a point to a center of a pixel (unit square) 
in the grid used by the Snap Rounding algorithm. Note that no conversion 
to an integer grid is done yet. Must have the syntax 
`void operator()(Point_2 p,FT pixel_size,FT &x,FT &y)` where \f$ p\f$ is the 
input point, `pixel_size` is the size of the pixel of the grid, 
and \f$ x\f$ and \f$ y\f$ are the \f$ x\f$ and \f$ y\f$-coordinates of the rounded point 
respectively. 
*/ 
typedef Hidden_type Snap_2; 

/*! 
Convert coordinates 
into an integer representation where one unit is equal to pixel size. 
For instance, if a point has the coordinates \f$ (3.7,5.3)\f$ and the pixel 
size is \f$ 0.5\f$, then the new point will have the coordinates of \f$ (7,10)\f$. 
Note, however, that the number type remains the same here, although 
integers are represented. 
Must have the syntax `Point_2 operator()(Point_2 p,NT pixel_size)` 
where \f$ p\f$ is the converted point and `pixel_size` is the size of the pixel 
of the grid. 
*/ 
typedef Hidden_type Integer_grid_point_2; 

/*! 
Returns the vertices of a polygon, 
which is the Minkowski sum of a segment and a square centered at the origin 
with edge size `pixel edge`. 
Must have the syntax 
`void operator()(std::list<Point_2>& vertices_list, Segment_2 s, 
NT unit_square)` 
where `vertices_list` is the list of the vertices of the Minkowski sum 
polygon, \f$ s\f$ is the input segment and `unit_square` is the edge size of 
the pixel. 
*/ 
typedef Hidden_type Minkowski_sum_with_pixel_2; 

/// @} 

/// \name Operations 
/// The following functions construct the required function objects
/// occasionally referred as functors listed above.
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

/// @}

}; /* end SnapRoundingTraits_2 */

