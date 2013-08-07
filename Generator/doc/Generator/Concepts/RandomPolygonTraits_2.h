/*!
\ingroup PkgGeneratorsConcepts
\cgalConcept

The concept `RandomPolygonTraits_2` describes the requirements for the traits 
class used by the function `random_polygon_2()`. 

\cgalHasModel \cgal kernels. 

\cgalHeading{Operations}

The following two member functions returning instances of the above predicate 
object types are required. 

*/

class RandomPolygonTraits_2 {
public:

/// \name Types 
/// @{

/*!
The coordinate type of the points of the polygon.
*/ 
typedef unspecified_type FT; 

/*!
The point type of the polygon. 
*/ 
typedef unspecified_type Point_2; 

/*!
Predicate object type that determines the orientation of three points. 
It must provide `Orientation operator()(Point_2 p, Point_2 q, 
Point_2 r)` that 
returns `LEFT_TURN`, if \f$ r\f$ lies to the left of the oriented 
line \f$ l\f$ defined by \f$ p\f$ and \f$ q\f$, returns `RIGHT_TURN` if \f$ r\f$ 
lies to the right of \f$ l\f$, and returns `COLLINEAR` if \f$ r\f$ lies 
on \f$ l\f$. 
*/ 
typedef unspecified_type Orientation_2; 

/*!
Binary predicate object type comparing `Point_2`s lexicographically. 
It must provide `bool operator()(Point_2 p, Point_2 q)` that 
returns `true` iff \f$ p <_{xy} q\f$. 
We have \f$ p<_{xy}q\f$, iff \f$ p_x < q_x\f$ or \f$ p_x = q_x\f$ and \f$ p_y < q_y\f$, 
where \f$ p_x\f$ and \f$ p_y\f$ denote the \f$ x\f$ and \f$ y\f$ coordinates of point \f$ p\f$, 
resp. 

*/ 
typedef unspecified_type Less_xy_2; 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Less_xy_2 less_xy_2_object(); 

/*!

*/ 
Orienation_2 orientation_2_object(); 

/// @}

}; /* end RandomPolygonTraits_2 */
