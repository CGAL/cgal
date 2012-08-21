/*!
\ingroup PkgPolygonPartitioning2Concepts
\cgalconcept

Function object that determines if a sequence of points represents a 
valid partition polygon or not, where "valid" can assume any of several 
meanings (e.g., convex or \f$ y\f$-monotone). 

\hasModel CGAL::Is_convex_2<Traits> 
\hasModel CGAL::Is_y_monotone_2<Traits> 

\sa `CGAL::approx_convex_partition_2` 
\sa `CGAL::convex_partition_is_valid_2` 
\sa `CGAL::greene_approx_convex_partition_2` 
\sa `CGAL::optimal_convex_partition_2` 
\sa `CGAL::partition_is_valid_2` 
\sa `CGAL::y_monotone_partition_2` 
\sa `CGAL::y_monotone_partition_is_valid_2` 

*/

class PolygonIsValid {
public:

/// \name Creation 
/// @{

/*! 

`Traits` is a model of the concept required by the function that checks 
for validity of the polygon. 

*/ 
PolygonIsValid(const Traits& t); 

/// @} 

/// \name Operations 
/// @{

/*! 

returns `true` iff the points of type `Traits::Point_2` 
in the range [`first`,`beyond`) define a valid polygon. 

*/ 
template<class InputIterator> 
bool operator()(InputIterator first, InputIterator beyond); 

/// @}

}; /* end PolygonIsValid */
