
/*!
\ingroup PkgInscribedAreasConcepts
\cgalConcept

\cgalAdvancedBegin
The concept `ExtremalPolygonTraits_2` provides the types and
operations needed to compute a maximal \f$ k\f$-gon that can be
inscribed into a given convex polygon.
\cgalAdvancedEnd

\note `ExtremalPolygonTraits_2::Less_xy_2` and 
`ExtremalPolygonTraits_2::Orientation_2` are used for (expensive) 
precondition checking only. Therefore, they need not to be 
specified, in case that precondition checking is disabled. 

\cgalHasModel `CGAL::Extremal_polygon_area_traits_2<K>`
\cgalHasModel `CGAL::Extremal_polygon_perimeter_traits_2<K>`

\sa `CGAL::maximum_area_inscribed_k_gon_2()`
\sa `CGAL::maximum_perimeter_inscribed_k_gon_2()`
\sa `CGAL::extremal_polygon_2()`

*/

class ExtremalPolygonTraits_2 {
public:

/// \name Types 
/// @{

/*!
model for `FieldNumberType`. 
*/ 
typedef unspecified_type FT; 

/*!
model for 
`Kernel::Point_2`. 
*/ 
typedef unspecified_type Point_2; 

/*!
model for 
`Kernel::Less_xy_2`. 
*/ 
typedef unspecified_type Less_xy_2; 

/*!
model for 
`Kernel::Orientation_2`. 
*/ 
typedef unspecified_type Orientation_2; 

/*!
AdaptableBinaryFunction class `op`: 
`Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ `FT`. 
Together with `init` this operation recursively defines the 
objective function to maximize. Let \f$ p\f$ and \f$ q\f$ be two vertices 
of a polygon \f$ P\f$ such that \f$ q\f$ precedes \f$ p\f$ in the oriented 
vertex chain of \f$ P\f$ starting with vertex \f$ root\f$. Then 
`op(p,q)` returns the value by which an arbitrary 
sub-polygon of \f$ P\f$ with vertices from \f$ [root,\, q]\f$ increases 
when \f$ p\f$ is added to it. E.g. in the maximum area case this is 
the area of the triangle \f$ (root,\, q,\, p)\f$. 
*/ 
typedef unspecified_type Operation; 

/// @} 

/// \name Operations 
/// @{

/*!
returns the minimal \f$ k\f$ for 
which a maximal \f$ k\f$-gon can be computed. (e.g. in the maximum 
area case this is three.) 
*/ 
int min_k() const; 

/*!
returns the value of the objective function for a 
polygon consisting of the two points `p` and `q`. (e.g. 
in the maximum area case this is `FT( 0)`.) 
*/ 
FT init( const Point_2& p, const Point_2& q) 
const; 

/*!
return `Operation` where `p` is the fixed \f$ root\f$ 
point. 
*/ 
Operation operation( const Point_2& p) 
const; 

/*!
writes the 
points of [`points_begin`, `points_end`) forming a 
`min_k()`-gon rooted at `points_begin[0]` of maximal 
value to o and returns the past-the-end iterator for that 
sequence (== `o + min_k()`). 
*/ 
template < class RandomAccessIterator, class 
OutputIterator > OutputIterator compute_min_k_gon( 
RandomAccessIterator points_begin, RandomAccessIterator 
points_end, FT& max_area, OutputIterator o) const; 

/*!

*/ 
Less_xy_2 less_xy_2_object(); 

/*!

*/ 
Orientation_2 orientation_2_object(); 

/// @}

}; /* end ExtremalPolygonTraits_2 */

