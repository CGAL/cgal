/*!
\ingroup PkgBarycentric_Coordinates_2Concepts
\cgalConcept

Requirements of the template parameter `Coordinate_2` for the class `CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2` with a set of two-dimensional generalized barycentric coordinates. 

\cgalHasModel `CGAL::Barycentric_coordinates::Wachspress_2`
\cgalHasModel `CGAL::Barycentric_coordinates::Mean_value_2`
\cgalHasModel `CGAL::Barycentric_coordinates::Discrete_harmonic_2`

*/

class BarycentricCoordinates_2 {

public:

/// \name Creation
/// @{

BarycentricCoordinates_2(const std::vector<Traits::Point_2> &vertices, const Traits &barycentric_traits);

/// @}

/// \name Functions
/// @{

/*!
	a function that computes generalized barycentric weights without normalization.
*/ 
std::pair<Iterator, bool> weights();

/*!
	a function that computes generalized barycentric coordinates on the bounded side of a polygon with an algorithm that is as precise as possible.
*/ 
std::pair<Iterator, bool> coordinates_on_bounded_side_precise();

/*!
	a function that computes generalized barycentric coordinates on the bounded side of a polygon with an algorithm that is as fast as possible.
*/ 
std::pair<Iterator, bool> coordinates_on_bounded_side_fast(); 

/*!
	a function that computes generalized barycentric coordinates on the unbounded side of a polygon with an algorithm that is as precise as possible.
*/ 
std::pair<Iterator, bool> coordinates_on_unbounded_side_precise();

/*!
	a function that computes generalized barycentric coordinates on the unbounded side of a polygon with an algorithm that is as fast as possible.
*/ 
std::pair<Iterator, bool> coordinates_on_unbounded_side_fast();

/*!
	a function that prints some information about computed coordinates.
*/ 
void print_coordinates_information();

/// @} 

}; /* end BarycentricCoordinates_2 */