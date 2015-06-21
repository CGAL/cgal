/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `FuzzyQueryItem` describes the requirements for fuzzy `d`-dimensional spatial objects. 

\cgalHasModel `CGAL::Fuzzy_sphere<Traits>`
\cgalHasModel `CGAL::Fuzzy_iso_box<Traits>` 

*/

class FuzzyQueryItem {
public:

/// \name Types 
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type D;

/*!
represents a `d`-dimensional point. 
*/ 
typedef unspecified_type Point_d; 

/*!
Number type. 
*/ 
typedef unspecified_type FT; 

/// @} 

/// \name Operations 
/// @{

/*!
Test whether the query item contains `p`. 
*/ 
bool contains(Point_d p) const; 

/*!
Test whether the inner approximation of the spatial object intersects a rectangle 
associated with a node of a tree. 
*/ 
bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const; 

/*!
Test whether the outer approximation of the spatial object encloses the rectangle 
associated with a node of a tree. 
*/ 
bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const; 

/// @}

}; /* end FuzzyQueryItem */
