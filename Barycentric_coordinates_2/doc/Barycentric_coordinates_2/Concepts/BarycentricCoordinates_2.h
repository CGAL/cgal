/*!
\ingroup PkgBarycentric_coordinates_2Concepts
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
	a function that computes generalized barycentric coordinates without normalization that are called generalized baycentric weights (as fast as possible algorithm is used)
*/ 
boost::optional<OutputIterator> weights(const Traits::Point_2 &query_point, OutputIterator &output);

/*!
	a function that computes generalized barycentric coordinates on the bounded side of a polygon with an algorithm that is as precise as possible
*/ 
boost::optional<OutputIterator> coordinates_on_bounded_side_precise(const Traits::Point_2 &query_point, OutputIterator &output);

/*!
	a function that computes generalized barycentric coordinates on the bounded side of a polygon with an algorithm that is as fast as possible
*/ 
boost::optional<OutputIterator> coordinates_on_bounded_side_fast(const Traits::Point_2 &query_point, OutputIterator &output); 

/*!
	a function that computes generalized barycentric coordinates on the unbounded side of a polygon with an algorithm that is as precise as possible
*/ 
boost::optional<OutputIterator> coordinates_on_unbounded_side_precise(const Traits::Point_2 &query_point, OutputIterator &output);

/*!
	a function that computes generalized barycentric coordinates on the unbounded side of a polygon with an algorithm that is as fast as possible
*/ 
boost::optional<OutputIterator> coordinates_on_unbounded_side_fast(const Traits::Point_2 &query_point, OutputIterator &output);

/// @} 

}; /* end BarycentricCoordinates_2 */
