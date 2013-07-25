
namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The class `Min_quadrilateral_default_traits_2` is a traits class for the 
functions `min_rectangle_2`, `min_parallelogram_2` and 
`min_strip_2` using a two-dimensional \cgal kernel. 

\tparam K must be a model for `Kernel`. 

\cgalModels `MinQuadrilateralTraits_2`

\sa `CGAL::min_rectangle_2()` 
\sa `CGAL::min_parallelogram_2()` 
\sa `CGAL::min_strip_2()` 

*/
template< typename K >
class Min_quadrilateral_default_traits_2 {
public:

/// \name Types 
/// @{

/*!
`Kernel::Point_2`. 
*/ 
typedef unspecified_type Point_2; 

/*!
`Kernel::Vector_2`. 
*/ 
typedef unspecified_type Vector_2; 

/*!
`Kernel::Direction_2`. 
*/ 
typedef unspecified_type Direction_2; 

/*!
`Kernel::Line_2`. 
*/ 
typedef unspecified_type Line_2; 

/*!
internal type. 
*/ 
typedef unspecified_type Rectangle_2; 

/*!
internal type. 
*/ 
typedef unspecified_type Parallelogram_2; 

/*!
internal type. 
*/ 
typedef unspecified_type Strip_2; 

/// @} 

/// \name Predicates 
/// @{

/*!
`Kernel::Equal_2`. 
*/ 
typedef unspecified_type Equal_2; 

/*!
`Kernel::Less_xy_2`. 
*/ 
typedef unspecified_type Less_xy_2; 

/*!
`Kernel::Less_yx_2`. 
*/ 
typedef unspecified_type Less_yx_2; 

/*!
`Kernel::Orientation_2`. 
*/ 
typedef unspecified_type Orientation_2; 

/*!
`Kernel::Has_on_negative_side_2`. 
*/ 
typedef unspecified_type Has_on_negative_side_2; 

/*!
`Kernel::Compare_angle_with_x_axis_2`. 
*/ 
typedef unspecified_type Compare_angle_with_x_axis_2; 

/*!
AdaptableBinaryFunction class 

`op`: `Rectangle_2` \f$ \times\f$ `Rectangle_2` 
\f$ \rightarrow\f$ `bool`. 
`op(r1,r2)` returns true, iff 
the area of \f$ r1\f$ is strictly less than the area of \f$ r2\f$. 
*/ 
typedef unspecified_type Area_less_rectangle_2; 

/*!
AdaptableBinaryFunction 
class 
`op`: `Parallelogram_2` \f$ \times\f$ 
`Parallelogram_2` \f$ \rightarrow\f$ `bool`. 

`op(p1,p2)` returns true, iff the area of \f$ p1\f$ is strictly 
less than the area of \f$ p2\f$. 
*/ 
typedef unspecified_type Area_less_parallelogram_2; 

/*!
AdaptableBinaryFunction class 

`op`: `Strip_2` \f$ \times\f$ `Strip_2` \f$ \rightarrow\f$ 
`bool`. 
`op(s1,s2)` returns true, iff the width of 
\f$ s1\f$ is strictly less than the width of \f$ s2\f$. 
*/ 
typedef unspecified_type Width_less_strip_2; 

/// @} 

/// \name Constructions 
/// @{

/*!
`Kernel::Construct_vector_2`. 
*/ 
typedef unspecified_type Construct_vector_2; 

/*!
AdaptableFunctor 

`op`: `Direction_2` \f$ \rightarrow\f$ `Vector_2`. 

`op(d)` returns a vector in direction `d`. 
*/ 
typedef unspecified_type Construct_vector_from_direction_2; 

/*!
`Kernel::Construct_perpendicular_vector_2`. 
*/ 
typedef unspecified_type Construct_perpendicular_vector_2; 

/*!
`Kernel::Construct_direction_2`. 
*/ 
typedef unspecified_type Construct_direction_2; 

/*!
`Kernel::Construct_opposite_direction_2`. 
*/ 
typedef unspecified_type Construct_opposite_direction_2; 

/*!
`Kernel::Construct_line_2`. 
*/ 
typedef unspecified_type Construct_line_2; 

/*!
Function class 
`op`: 
`Point_2` \f$ \times\f$ `Direction_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`Rectangle_2`. 
If the points 
`p1`,`p2`,`p3`,`p4` form the boundary of a 
convex polygon (oriented counterclockwise), 
`op(p1,d,p2,p3,p4)` returns the rectangle with one of the 
points on each side and one sides parallel to `d`. 
*/ 
typedef unspecified_type Construct_rectangle_2; 

/*!
Function class 

`op`: `Point_2` \f$ \times\f$ `Direction_2` \f$ \times\f$ 
`Point_2` \f$ \times\f$ `Direction_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `Point_2` \f$ \rightarrow\f$ `Rectangle_2`. 
If the 
points `p1`,`p2`,`p3`,`p4` form the 
boundary of a convex polygon (oriented counterclockwise), 
`op(p1,d1,p2,d2,p3,p4)` returns the parallelogram with one 
of the points on each side and one side parallel to each of 
`d1` and `d2`. 
*/ 
typedef unspecified_type Construct_parallelogram_2; 

/*!
Function class 
`op`: 
`Point_2` \f$ \times\f$ `Direction_2` \f$ \times\f$ `Point_2` 
\f$ \rightarrow\f$ `Strip_2`. 
`op(p1,d,p2)` returns the 
strip bounded by the lines through `p1` resp. `p2` with 
direction `d`. 
*/ 
typedef unspecified_type Construct_strip_2; 

/// @} 

/// \name Operations 
/// Additionally, for each of the predicate and construction functor
/// types listed above, there is a member function that requires no
/// arguments and returns an instance of that functor type. The name
/// of the member function is the uncapitalized name of the type
/// returned with the suffix `_object` appended. For example, for the
/// functor type `Construct_vector_2` the following member function
/// exists:
/// @{

/*!
copies the four vertices of `r` in 
counterclockwise order to `o`. 
*/ 
template < class OutputIterator > OutputIterator 
copy_rectangle_vertices_2(const Rectangle_2& r, OutputIterator 
o) const; 

/*!
copies the four vertices of `p` in 
counterclockwise order to `o`. 
*/ 
template < class OutputIterator > OutputIterator 
copy_parallelogram_vertices_2(const Parallelogram_2& p, 
OutputIterator o) const; 

/*!
copies the two lines bounding `s` to `o`. 
*/ 
template < class OutputIterator > OutputIterator 
copy_strip_lines_2(const Strip_2& s, OutputIterator o) 
const; 

/*!

*/ 
Construct_vector_2 construct_vector_2_object() 
const ; 

/// @}

}; /* end Min_quadrilateral_default_traits_2 */
} /* end namespace CGAL */
