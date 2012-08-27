/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalconcept

Requirements of an orthogonal distance class supporting incremental distance updates. 
To optimize distance computations transformed distances are used. 
E.g., for an Euclidean distance the transformed distance is the squared Euclidean distance. 

\refines ::GeneralDistance 

\hasModel `CGAL::Euclidean_distance<Traits>`
\hasModel `CGAL::Weighted_Minkowski_distance<Traits>`

*/

class OrthogonalDistance {
public:

/// \name Types 
/// @{

/*! 
Number type. 
*/ 
typedef Hidden_type FT; 

/*! 
Point type. 
*/ 
typedef Hidden_type Point_d; 

/*! 
Query item type. 
*/ 
typedef Hidden_type Query_item; 

/// @} 

/// \name Creation 
/// @{

/*! 
Constructor implementing distance for \f$ d\f$-dimensional points. 
*/ 
OrthogonalDistance(int d); 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the transformed distance between `q` and `r`. 
*/ 
FT transformed_distance(Query_item q, Point_d r) const; 

/*! 
Returns the transformed distance between `q` and 
the point on the boundary of `r` closest to `q`. 
*/ 
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT> r) const; 

/*! 
Returns the transformed distance between `q` and 
the point on the boundary of `r` farthest to `q`. 
*/ 
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT> r) const; 

/*! 
Returns the transformed distance. 
*/ 
FT transformed_distance(FT d) const; 

/*! 
Returns the inverse of the transformed distance. 
*/ 
FT inverse_of_transformed_distance(FT d) const; 

/*! 
Updates 
`dist` incrementally and returns the updated distance. 
*/ 
FT new_distance(FT dist, FT old_off, FT new_off, int cutting_dimension) const; 

/// @}

}; /* end OrthogonalDistance */
