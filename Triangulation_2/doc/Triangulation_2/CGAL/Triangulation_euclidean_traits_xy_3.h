
namespace CGAL {

/*!
\ingroup PkgTriangulation2TraitsClasses

The functionality of this class has been generalized to other packages than 2D triangulations. 
The more general class `Projection_traits_xy_3` can be found in the 2D and 3D Linear Geometric Kernel. 

\deprecated The class `Triangulation_euclidean_traits_xy_3` is a
geometric traits class which allows to triangulate a terrain. This
traits class is designed to build a two dimensional triangulation
embedded in 3D space, i.e.\ a triangulated surface, such that its on
the \f$ xy\f$ plane is a Delaunay triangulation.  This is a usual
construction for GIS terrains.  Instead of really projecting the 3D
points and maintaining a mapping between each point and its projection
(which costs space and is error prone) the class
`Triangulation_euclidean_traits_xy_3` supplies geometric predicates
that ignore the `z`-coordinate of the points.

The class is a model of the concept `DelaunayTriangulationTraits_2` 
except that it does not provide the type and constructors 
required to build the dual Voronoi diagram. The class is also a model 
of the concept `ConstrainedTriangulationTraits_2`. 

\cgalHeading{Parameters}

The template parameter `K` has to 
be instantiated by a model of the `Kernel` concept. 
`Triangulation_euclidean_traits_xy_3` uses types 
and predicates defined in `K`. 

\sa `TriangulationTraits_2` 
\sa `DelaunayTriangulationTraits_2` 
\sa `CGAL::Triangulation_2<Traits,Tds>` 
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>` 

\cgal provides also predefined geometric traits class
`Triangulation_euclidean_traits_yz_3<K>` and
`Triangulation_euclidean_traits_xz_3<K>` to deal with projections on
the `xz`- or the `yz`-plane, respectively.

\sa `CGAL/Triangulation_euclidean_traits_xz_3.h`
\sa `CGAL/Triangulation_euclidean_traits_yz_3.h`

*/
template< typename K >
class Triangulation_euclidean_traits_xy_3 {
public:

/// \name Types 
/// The following predicates and constructor types are provided
/// @{

/*!

*/ 
typedef Point_3<K> Point_2; 

/*!

*/ 
typedef Segment_3<K> Segment_2; 

/*!

*/ 
typedef Triangle_3<K> Triangle_2; 

/*!

*/ 
typedef Line_3<K> Line_2; 

/*!
A constructor object for 
`Segment_2`. Provides : 

`Segment_2 operator()(Point_2 p,Point_2 q)`, 

which constructs a segment from two points. 
*/ 
typedef unspecified_type Construct_segment_2; 

/*!
A constructor object for 
`Triangle_2`. Provides : 

`Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r )`, 

which constructs a triangle from three points. 
*/ 
typedef unspecified_type Construct_triangle_2; 

/*!
A constructor object for 
`Line_2`. Provides : 

`Line_2 operator()(Point_2 p,Point_2 q)`, 

which constructs a line from two points. 
*/ 
typedef unspecified_type Construct_line_2; 

/*!
A construction object. 
Provides the operator : 

`RT operator()(Line_2 l, Point_2 p);` 
which returns the squared distance between the projection of `p` 
and the projection of `l`. 
*/ 
typedef unspecified_type Compute_squared_distance_2; 

/*!
A construction object. 
Provides the operator : 

`Object_2 operator()(Segment_2 s1, Segment_2 s2);` 
which returns the intersection of the projection of `s1` 
and the projection of `s2` embedded in `3D`. If the intersection 
is a segment, the `z`-coordinates of its extremities is `0`. 
If the intersection is a point `p`, let `p1` and `p2` be the points on `s1` 
and `s2` respectively, such that their projections are `p`. The point returned is the 
middle of the segment `p1``p2`. 
\pre The projection of `s1` and the projection of `s2` are non-degenerate `2D` segments. 

*/ 
typedef unspecified_type Intersect_2; 

/*!
Predicate object. Provides 
the operator : 

`Comparison_result operator()(Point_2 p, Point_2 q)` 

which returns 
`SMALLER, EQUAL` or `LARGER` 
according to the 
\f$ x\f$-ordering of points `p` and `q`. 
*/ 
typedef unspecified_type Compare_x_2; 

/*!
Predicate object. Provides 
the operator : 

`Comparison_result operator()(Point_2 p, Point_2 q)` 

which returns 
(`SMALLER, EQUAL` or `LARGER`) 
according to the 
\f$ y\f$-ordering of points `p` and `q`. 
*/ 
typedef unspecified_type Compare_y_2; 

/*!
Predicate object. Provides 
the operator : 

`Orientation operator()(Point_2 p, Point_2 q, Point_2 r)` 

which returns 
`LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR` 
according to the position of the projection of \f$ r\f$ 
with respect to the projection of the 
oriented line `pq`. 
*/ 
typedef unspecified_type Orientation_2; 

/*!
Predicate object. 
Provides the operator : 
`Oriented_side operator()(Point_2 p, Point_2 q, Point_2 r, Point_2 s)` 
which takes four points \f$ p, q, r, s\f$ as arguments and returns 
`ON_POSITIVE_SIDE`, `ON_NEGATIVE_SIDE` or, 
`ON_ORIENTED_BOUNDARY` according to the position of 
the projection of point`s` 
with respect to the oriented circle through the projections of \f$ p,q\f$ 
and \f$ r\f$. 
*/ 
typedef unspecified_type Side_of_oriented_circle_2; 

/*!
Predicate object. Provides 
the operator : 

`Comparison_result operator()(Point_2 p, Point_2 q, Point_2 r)` 
which returns `SMALLER`, `EQUAL` or `LARGER` 
according to the distance between the projection of p and the projection of q being smaller, equal or larger 
than the distance between the projection of p and the projection of r. 
*/ 
typedef unspecified_type Compare_distance_2; 

/// @} 

/// \name Creation 
/// @{

/*!

default constructor. 
*/ 
Triangulation_euclidean_traits_xy_3(); 

/*!
Copy constructor. 
*/ 
Triangulation_euclidean_traits_xy_3( 
Triangulation_euclidean_traits_xy_3 tr); 

/*!
Assignment operator. 
*/ 
Triangulation_euclidean_traits_xy_3 operator= 
(Triangulation_euclidean_traits_xy_3 tr); 

/// @} 

/// \name Access to predicate objects 
/// The following access functions are provided 
/// @{

/*!

*/ 
Construct_segment_2 construct_segment_2_object(); 

/*!

*/ 
Construct_triangle_2 construct_triangle_2_object(); 

/*!

*/ 
Construct_line_2 construct_line_2_object(); 

/*!

*/ 
Comparison_x_2 compare_x_2_object(); 

/*!

*/ 
Comparison_y_2 compare_y_2_object(); 

/*!

*/ 
Orientation_2 orientation_2_object(); 

/*!

*/ 
Side_of_oriented_circle_2 
side_of_oriented_circle_2_object(); 

/*!

*/ 
Compare_distance_2 
compare_distance_2_object(); 

/*!

*/ 
Intersect_2 intersect_2_object(); 

/*!

*/ 
Compute_squared_distance_2 compute_squared_distance_2_object(); 

/// @}

}; /* end Triangulation_euclidean_traits_xy_3 */
} /* end namespace CGAL */
