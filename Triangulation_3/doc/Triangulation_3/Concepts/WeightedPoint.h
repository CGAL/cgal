
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

\cgalHasModel `CGAL::Weighted_point`

\sa `CGAL::Regular_triangulation_euclidean_traits_3` 
\sa `CGAL::Regular_triangulation_3`

The concept `WeightedPoint` is needed by
`CGAL::Regular_triangulation_euclidean_traits_3`.
It must fulfill the following requirements:

*/

class WeightedPoint {
public:

/// \name Types 
/// @{

/*!
The point type 
*/ 
typedef unspecified_type Point; 

/*!
The weight type 
*/ 
typedef unspecified_type Weight; 

/*!
The ring type 
*/ 
typedef Point::RT RT; 

/// @} 

/// \name Creation 
/// @{

/*!

*/ 
WeightedPoint(const Point &p=Point(), const Weight &w 
= Weight(0)); 

/// @}

/// \name Access Functions
/// @{

/*!

*/ 
Point point() const; 

/*!

*/ 
Weight weight() const; 

/// @}

}; /* end WeightedPoint */

