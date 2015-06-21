/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

Requirements of a distance class defining a distance between a query item 
denoting a spatial object and a point. 
To optimize distance computations transformed distances are used, 
e.g., for a Euclidean distance the transformed distance is the squared 
Euclidean distance. 

\cgalHasModel `CGAL::Manhattan_distance_iso_box_point<Traits>` 
\cgalHasModel `CGAL::Euclidean_distance_sphere_point<Traits>`

*/

class GeneralDistance {
public:

/// \name Types 
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type D;

/*!
Number type. 
*/ 
typedef unspecified_type FT; 

/*!
Point type. 
*/ 
typedef unspecified_type Point_d; 

/*!
Query item type. 
*/ 
typedef unspecified_type Query_item; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the transformed distance between `q` and `r`. 
*/ 
FT transformed_distance(Query_item q, Point_d r); 

/*!
Returns the transformed distance between `q` and 
the point on the boundary of `r` closest to `q`. 
*/ 
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the transformed distance between `q` and 
the point on the boundary of `r` furthest to `q`. 
*/ 
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the transformed distance. 
*/ 
FT transformed_distance(FT d) const; 

/*!
Returns the inverse of the transformed distance. 
*/ 
FT inverse_of_transformed_distance(FT d) const; 

/// @}

}; /* end GeneralDistance */
