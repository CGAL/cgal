/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

All convex hull and extreme point algorithms provided in \cgal are 
parameterized with a traits class `Traits`, which defines the 
primitives (objects and predicates) that the cone based spanner algorithms use. 

\cgalHasModel any model of a \cgal %kernel.

*/

class ConeBasedSpannerTraits {
public:

/// \name Types 
/// @{

/*!
The point type. 
*/ 
typedef unspecified_type Point_2; 


/*!
Predicate object type that must provide 
`bool operator()(Point_2 p, Point_2 q, 
Point_2 r,Point_2 s)`, which returns `true` iff 
the signed distance from \f$ r\f$ to the line \f$ l_{pq}\f$ through \f$ p\f$ and \f$ q\f$ 
is smaller than the distance from \f$ s\f$ to \f$ l_{pq}\f$. It is used to 
compute the point right of a line with maximum unsigned distance to 
the line. The predicate must provide a total order compatible 
with convexity, <I>i.e.</I>, for any line segment \f$ s\f$ one of the 
endpoints 
of \f$ s\f$ is the smallest point among the points on \f$ s\f$, with respect to 
the order given by `Less_signed_distance_to_line_2`. 
*/ 
typedef unspecified_type Less_signed_distance_to_line_2; 


/// @} 



/// \name Operations 
/// The following member functions to create instances of the above predicate 
/// object types must exist. 
/// @{


/*!

*/ 
Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object(); 


/// @}

}; /* end ConeBasedSpannerTraits */
