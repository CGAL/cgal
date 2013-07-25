
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `GpsTraitsGeneralPolygonWithHoles_2` 

A model of this concept represents a general polygon with holes. The 
concept requires the ability to access the general polygon that 
represents the outer boundary and the general polygons that represent 
the holes. 

\cgalHasModel `CGAL::General_polygon_with_holes_2<General_polygon>`
\cgalHasModel `CGAL::Polygon_with_holes_2<Kernel,Container>`
\cgalHasModel `CGAL::Gps_circle_segment_traits_2<Kernel>::%Polygon_with_holes_2`
\cgalHasModel `CGAL::Gps_traits_2<ArrTraits,GeneralPolygon>::%Polygon_with_holes_2` 

*/

class GeneralPolygonWithHoles_2 {
public:

/// \name Types 
/// @{

/*!
the general-polygon type used to 
represent the outer boundary and each hole. Must model the `GeneralPolygon_2` concept. 
*/ 
typedef unspecified_type General_polygon_2; 

/*!
a bidirectional iterator 
over the polygonal holes. Its value type is 
`General_polygon_2`. 
*/ 
typedef unspecified_type Hole_const_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
GeneralPolygonWithHoles_2(); 

/*!
copy constructor. 
*/ 
GeneralPolygonWithHoles_2(GeneralPolygonWithHoles_2 other); 

/*!
assignment operator. 
*/ 
GeneralPolygonWithHoles_2 operator=(other); 

/*!
constructs a general polygon with holes that has no holes using a given general polygon `outer` as the outer boundary. 
*/ 
GeneralPolygonWithHoles_2(General_polygon_2 & outer); 

/*!
constructs a general polygon with holes using a given general polygon 
`outer` as the outer boundary and a given range of holes. If `outer` 
is an empty general polygon, then an unbounded polygon with holes will be 
created. The holes must be contained inside the outer boundary, and the 
polygons representing the holes must be simple and pairwise disjoint, except 
perhaps at the vertices. 
*/ 
template <class InputIterator> 
GeneralPolygonWithHoles_2(General_polygon_2 & outer, 
InputIterator begin, InputIterator end); 

/// @} 

/// \name Predicates 
/// @{

/*!
returns `true` if the outer boundary is empty, and `false` 
otherwise. 
*/ 
bool is_unbounded(); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the general polygon that represents the outer boundary. Note that this polygon is not necessarily a valid (simple) general polygon because it may be relatively simple. 
*/ 
const General_polygon_2 & outer_boundary() const; 

/*!
returns the begin iterator of the holes. 
*/ 
Hole_const_iterator holes_begin() const; 

/*!
returns the past-the-end iterator of the holes. 
*/ 
Hole_const_iterator holes_end() const; 

/// @}

}; /* end GeneralPolygonWithHoles_2 */

