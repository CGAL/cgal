/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalconcept

The concept `FuzzyQueryItem` describes the requirements for fuzzy \f$ d\f$-dimensional spatial objects. 

\hasModel CGAL::Fuzzy_sphere<Traits>
\hasModel CGAL::Fuzzy_iso_box<Traits> 

*/

class FuzzyQueryItem {
public:

/// \name Types 
/// @{

/*! 
represents a \f$ d\f$-dimensional point. 
*/ 
typedef Hidden_type Point_d; 

/*! 
Number type. 
*/ 
typedef Hidden_type FT; 

/// @} 

/// \name Operations 
/// @{

/*! 
test whether \f$ q\f$ contains \f$ p\f$. 
*/ 
bool contains(Point_d p) const; 

/*! 
test whether the inner approximation of the spatial object intersects a rectangle 
associated with a node of a tree. 
*/ 
bool inner_range_intersects(const Kd_tree_rectangle<FT>& rectangle) const; 

/*! 
test whether the outer approximation of the spatial object encloses the rectangle 
associated with a node of a tree. 
*/ 
bool outer_range_contains(const Kd_tree_rectangle<FT>& rectangle) const; 

/// @}

}; /* end FuzzyQueryItem */
