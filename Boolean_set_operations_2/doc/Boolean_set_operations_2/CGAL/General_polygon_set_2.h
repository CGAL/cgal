namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

An object of the `General_polygon_set_2` class-template represents a 
point set in the plane bounded by \f$ x\f$ monotone curves. Points in the set 
lie on the boundary or on the positive side of the curves. This class 
template provides methods to apply regularized Boolean set-operations and 
few other utility methods. An `Arrangement_2` data structure is 
internally used to represent the point set. The arrangement is 
represented as a doubly-connected edge-list (<span class="textsc">Dcel</span>). 

The `Traits` template-parameter should be instantiated with a 
model of the concept `GeneralPolygonSetTraits_2`. The traits class 
defines the types of points, \f$ x\f$-monotone curves, general polygons, 
and general polygons with holes, that is 
`Traits::Point_2`, 
`Traits::X_monotone_curve_2`, `Traits::Polygon_2`, and 
`Traits::Polygon_with_holes_2`, respectively. 
`Traits::Point_2` must 
be the type of the endpoints of 
`Traits::X_monotone_curve_2`, and 
`Traits::X_monotone_curve_2` must be the type of 
the curves that comprise the boundaries of the general polygons. The traits 
class supports geometric operations on the types above. We sometimes use 
the term <I>polygon</I> instead of general polygon for simplicity hereafter. 

The template parameter `Dcel` should be instantiated with a 
model of the concept `GeneralPolygonSetDcel`. It is instantiated 
by default with the type `Gps_default_dcel<Traits>`. You can override 
this default, with a different <span class="textsc">Dcel</span> class, typically an extension 
of the `Gps_default_dcel` class template. Overriding the default is 
necessary only if you intend to obtain the underlying internal arrangement 
and process it further. 

The input and output of the Boolean set-operations methods consist of one 
or more general polygons, some of which may have holes. In particular, 
these methods operate on pairs of objects of type `General_polygon_set_2`, or 
directly on objects of type `Traits::Polygon_2` or 
`Traits::Polygon_with_holes_2`. An object of type 
`Traits::Polygon_2` is a valid operand, only if it 
is simple and its boundary is oriented counterclockwise. An object of type 
`Traits::Polygon_with_holes_2` is valid, only if 
its outer boundary bounds a relatively simple general polygon oriented 
counterclockwise, and each one of its holes is a simple polygon oriented 
clockwise. The holes are pairwise disjoint, except perhaps at the vertices, 
and are contained in the polygon bounded by the outer boundary. The outer 
boundary and the holes are also pairwise disjoint, except perhaps at the 
vertices. 

\sa `Arrangement_2` 
\sa `ArrangementXMonotoneTraits_2` 
\sa `Nef_polyhedron_2` 

*/
template< typename Traits, typename Dcel >
class General_polygon_set_2 {
public:

/// \name Types 
/// @{

/*!
the traits class in use. 
*/ 
typedef unspecified_type Traits_2; 

/*!
the general polygon type. Must model the `GpsTraitsGeneralPolygon_2` concept. 
*/ 
typedef unspecified_type Polygon_2; 

/*!
the general polygon with holes type. Must model the `GpsTraitsGeneralPolygonWithHoles_2` concept. 
*/ 
typedef unspecified_type Polygon_with_holes_2; 

/*!
number of polygons with holes size type. 
*/ 
typedef unspecified_type Size; 

/*!
the arrangement type used internally. 
*/ 
typedef unspecified_type Arrangement_2; 

/// @} 

/// \name Creation 
/// @{

/*!
constructs an empty set of polygons `gps` represented by an empty arrangement. 
*/ 
General_polygon_set_2<Traits>(); 

/*!
copy constructor. 
*/ 
General_polygon_set_2<Traits>(const Self & other); 

/*!
constructs an empty set of polygons  `gps` that uses the given 
`traits` instance for performing the geometric operations. 
*/ 
General_polygon_set_2<Traits>(Traits & traits); 

/*!
constructs a set of polygons  `gps` that consists of the single polygon `pgn`. 
*/ 
General_polygon_set_2<Traits>(const Polygon_2 & pgn); 

/*!
constructs a set of polygons  `gps` that consists of the single polygon with 
holes `pgn_with_holes`. 
*/ 
General_polygon_set_2<Traits>(const Polygon_with_holes_2 & pgn_with_holes); 

/// @} 

/// \name Access Functions 
/// @{

/*!
obtains the polygons with holes represented by `gps`. 
*/ 
template <class OutputIterator> 
OutputIterator polygons_with_holes(OutputIterator out); 

/*!
returns the total number of general polygons represented by `gps`. 
*/ 
Size number_of_polygons_with_holes() const; 

/*!
returns `true` if `gps` represents an empty set. 
*/ 
bool is_empty() const; 

/*!
returns `true` if `gps` represents the entire plane. 
*/ 
bool is_plane() const; 

/*!
obtains an instance of the traits. If the traits was passed as a 
parameter to the constructor of `gps`, it is returned. Otherwise, a 
newly created instance is returned. 
*/ 
Traits & traits() const; 

/*!
obtains the arrangement data structure that internally represents the 
general-polygon set. 
*/ 
const Arrangement_2 & arrangement() const; 

/// @} 

/// \name Modifiers 
/// @{

/*!
clears `gps`. 
*/ 
void clear(); 

/*!
inserts `pgn` into `gps`. 
\pre `pgn` and the point set represented by `gps` are disjoint. 
This precondition enables the use of very efficient insertion methods. 
Use the respective overloaded method that inserts a polygon of type 
`Polygon_with_holes_2`, if only a relaxed 
precondition can be guaranteed. If even the relaxed precondition cannot 
be guaranteed, use the `join` method. 
*/ 
void insert(Polygon_2 & pgn); 

/*!
inserts `pgn_with_holes` into `gps`. 
\pre `pgn_with_holes` does not intersect with the point set represented by `gps`, except maybe at the vertices. 
If this relaxed precondition cannot be guaranteed, use the `join` 
method. 
*/ 
void insert(Polygon_with_holes_2 & pgn_with_holes); 

/*!
inserts the range of polygons (or polygons with holes) into `gps`. (The 
value type of the input iterator is used to distinguish between the two.) 
\pre If the given range contains objects of type `Polygon_with_holes_2`, then these polygons and the point set represented by `gps` are pairwise disjoint, except maybe at the vertices. If the given range contains objects of type `Polygon_2`, then these polygons and the point set represented by `gps` are pairwise disjoint without any exceptions. 
*/ 
template <class InputIterator> 
void insert(InputIterator begin, InputIterator end); 

/*!
inserts the two ranges of polygons and polygons with holes into `gps`. 
\pre All polygons in the first range, all polygon with holes in the second range, and the point set represented by `gps` are pairwise disjoint, except maybe at the vertices 
*/ 
template <class InputIterator1, class InputIterator2> 
void insert(InputIterator1 pgn_begin, InputIterator1 pgn_end, 
InputIterator2 pgn_with_holes_begin, 
InputIterator2 pgn_with_holes_end); 

/*!
computes the complement of `gps`. 
*/ 
void complement(); 

/*!
computes the complement of `other`. 
`gps` is overridden by the result. 
*/ 
void complement(const Polygon_set_2 & other); 

/// @} 

/// \name Univariate Operations 
/// In the following univariate and bivariate methods the result is placed in `gps` after it is cleared.
/// @{

/*!
computes the intersection of `gps` and `other`. 
*/ 
void intersection(const General_polygon_set_2 & other); 

/*!
computes the intersection of `gps` and `pgn`. 
*/ 
void intersection(const Polygon_2 & pgn); 

/*!
computes the intersection of `gps` and `pgn`. 
*/ 
void intersection(const Polygon_with_holes_2 & pgn); 

/*!
computes the intersection of a collection of point sets. The collection 
consists of the polygons (or polygons with holes) in the given range and 
the point set represented by `gps`. (The value type of the input iterator 
is used to distinguish between the two options.) 
*/ 
template <class InputIterator> 
void intersection(InputIterator begin, InputIterator end); 

/*!
computes the intersection of a collection of point sets. The collection 
consists of the polygons and polygons with holes in the given two ranges 
and the point set represented by `gps`. 
*/ 
template <class InputIterator1, class InputIterator2> 
void intersection(InputIterator1 pgn_begin, InputIterator1 pgn_end, 
InputIterator2 pgn_with_holes_begin, 
InputIterator2 pgn_with_holes_end); 

/*!
computes the union of `gps` and `other`. 
*/ 
void join(const General_polygon_set_2 & other); 

/*!
computes the union of `gps` and `pgn`. 
*/ 
void join(const Polygon_2 & pgn); 

/*!
computes the union of `gps` and `pgn`. 
*/ 
void join(const Polygon_with_holes_2 & pgn); 

/*!
computes the union of the polygons (or polygons with holes) in the 
given range and the point set represented by `gps`. (The value type 
of the input iterator is used to distinguish between the two options.) 
*/ 
template <class InputIterator> 
void join(InputIterator begin, InputIterator end); 

/*!
computes the union of the polygons and polygons with holes in the 
given two ranges and the point set represented by `gps`. 
*/ 
template <class InputIterator1, class InputIterator2> 
void join(InputIterator1 pgn_begin, InputIterator1 pgn_end, 
InputIterator2 pgn_with_holes_begin, 
InputIterator2 pgn_with_holes_end); 

/*!
computes the difference between `gps` and `other`. 
*/ 
void difference(const General_polygon_set_2 & other); 

/*!
computes the difference between `gps` and `pgn`. 
*/ 
void difference(const Polygon_2 & pgn); 

/*!
computes the difference between `gps` and `pgn`. 
*/ 
void difference(const Polygon_with_holes_2 & pgn); 

/*!
computes the symmetric difference between `gps` and `other`. 
*/ 
void symmetric_difference(const General_polygon_set_2 & other); 

/*!
computes the symmetric difference between `gps` and `pgn`. 
*/ 
void symmetric_difference(const Polygon_2 & pgn); 

/*!
computes the symmetric difference between `gps` and `pgn`. 
*/ 
void symmetric_difference(const Polygon_with_holes_2 & pgn); 

/*!
computes the symmetric difference (xor) of a collection of point sets. 
The collection consists of the polygons (or polygons with holes) in the 
given range and the point set represented by `gps`. (The value type of 
the input iterator is used to distinguish between the two options.) 
*/ 
template <class InputIterator> 
void symmetric_difference(InputIterator begin, InputIterator end); 

/*!
computes the symmetric difference (xor) of a collection of point sets. 
The collection consists of the polygons and polygons with holes in the 
given two ranges and the point set represented by `gps`. 
*/ 
template <class InputIterator1, class InputIterator2> 
void symmetric_difference(InputIterator1 pgn_begin, InputIterator1 pgn_end, 
InputIterator2 pgn_with_holes_begin, 
InputIterator2 pgn_with_holes_end); 

/// @} 

/// \name Bivariate Operations 
/// The following bivariate function replace `gps` with the result.
/// @{

/*!
computes the intersection of `gps1` and `gps2`. 
*/ 
void intersection(const General_polygon_set_2 & gps1, 
const General_polygon_set_2 & gps2); 

/*!
computes the union of `gps1` and `gps2`. 
*/ 
void join(const General_polygon_set_2 & gps1, 
const General_polygon_set_2 & gps2); 

/*!
computes the difference between `gps1` and `gps2`. 
*/ 
void difference(const General_polygon_set_2 & gps1, 
const General_polygon_set_2 & gps2); 

/*!
computes the symmetric difference between `gps1` and `gps2`. 
*/ 
void symmetric_difference(const General_polygon_set_2 & gps1, 
const General_polygon_set_2 & gps2); 

/// @} 

/// \name Query Functions 
/// @{

/*!
returns `true` if `gps` and `other` intersect in their 
interior, and `false` otherwise. 
*/ 
bool do_intersect(const General_polygon_set_2 & other); 

/*!
returns `true` if `gps` and `pgn` intersect in their 
interior, and `false` otherwise. 
*/ 
bool do_intersect(const Polygon_2 & pgn); 

/*!
returns `true` if `gps` and `pgn` intersect in their 
interior, and `false` otherwise. 
*/ 
bool do_intersect(const Polygon_with_holes_2 & pgn); 

/*!
returns `true` if the interior of the point sets in a collection 
intersect, and `false` otherwise. The collection consists of the 
polygons (or polygons with holes) in the given range and the point set 
represented by `gps`. (The value type of the input iterator is used to 
distinguish between the two options.) 
*/ 
template <class InputIterator> 
void do_intersect(InputIterator begin, InputIterator end); 

/*!
returns `true` if the interior of the point sets in a collection 
intersect, and `false` otherwise. The collection consists of the 
polygons and polygons with holes in the given two ranges and the point 
set represented by `gps`. 
*/ 
template <class InputIterator1, class InputIterator2> 
void do_intersect(InputIterator1 pgn_begin, InputIterator1 pgn_end, 
InputIterator2 pgn_with_holes_begin, 
InputIterator2 pgn_with_holes_end); 

/*!
obtains a polygon with holes that contains the query point `p`, 
if exists, through `pgn`, and returns `true`. 
Otherwise, returns `false`. 
*/ 
bool locate(const Point_2 & p, Polygon_with_holes_2 & pgn); 

/*!
returns either the constant `ON_ORIENTED_BOUNDARY`, 
`ON_POSITIVE_SIDE`, or `ON_NEGATIVE_SIDE`, iff `p` lies on 
the boundary, properly on the positive side, or properly on the negative 
side of `gps` respectively. 
*/ 
Oriented_side oriented_side(const Point_2 & q); 

/*!
returns either the constant `ON_NEGATIVE_SIDE`, 
`ON_ORIENTED_BOUNDARY`, or `ON_POSITIVE_SIDE`, iff 
`other` and `gps` are completely disjoint, in contact, or 
intersect in their interior, respectively. 
*/ 
Oriented_side oriented_side(const General_polygon_set_2 & other); 

/*!
returns either the constant `ON_NEGATIVE_SIDE`, 
`ON_ORIENTED_BOUNDARY`, or `ON_POSITIVE_SIDE`, iff 
`pgn` and `gps` are completely disjoint, in contact, or 
intersect in their interior, respectively. 
*/ 
Oriented_side oriented_side(const Polygon_2 & pgn); 

/*!
returns either the constant `ON_NEGATIVE_SIDE`, 
`ON_ORIENTED_BOUNDARY`, or `ON_POSITIVE_SIDE`, iff 
`pgn` and `gps` are completely disjoint, in contact, or 
intersect in their interior, respectively. 
*/ 
Oriented_side oriented_side(const Polygon_with_holes_2 & pgn); 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
returns `true` if `gps` represents a valid point set. 
*/ 
bool is_valid() const; 

/// @}

}; /* end General_polygon_set_2 */
} /* end namespace CGAL */
