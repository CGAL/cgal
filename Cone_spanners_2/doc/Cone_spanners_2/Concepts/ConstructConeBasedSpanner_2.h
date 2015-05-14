
/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The concept `ConeSpanners_2` describes the set of requirements to be fulfilled by any class that constructs a cone based spanner, such as a Theta graph and a Yao graph, etc. 

\cgalHasModel `CGAL::Theta_graph_2`
\cgalHasModel `CGAL::Yao_graph_2`

*/

class ConeSpanners_2 {
public:

/// \name Types 
/// @{

/*!
Any model of this concept must use `boost::adjacency_list` to store the constructed cone based spanner. Of the seven
template parameters of `boost::adjacency_list`, the fourth parameter 'VertexProperties' must be Point_2 from \cgal, 
and other parameters can be decided freely. This particular type is named 'Graph'.
*/
typedef adjacency_list<OutEdgeList, VertexList, Directed,
               Point_2, EdgeProperties,
               GraphProperties, EdgeList> Graph;

/// @}


/// \name Creation 
/// @{

/** @brief constructor for a ConeSpanners_2 object.
*
* @tparam PointInputIterator   Must be Point_2 of a certain type of \cgal Kernel.
* @param k     Number of cones to divide space into
* @param start An iterator pointing to the first point (vertex) in the graph.
*              (default: nullptr)
* @param end   An iterator pointing to the place that passes the last point.  (default: nullptr)
* @param ray0  The direction of the first ray. This allows the first ray to be at an arbitary 
*              direction.  (default: positive x-axis) 
*/
template <typename PointInputIterator>  
ConeSpanners_2(const unsigned int k,
			const PointInputIterator& start=nullptr, 
			const PointInputIterator& end=nullptr,
			const Direction_2& ray0 = Direction_2(1,0)	);

/// @} 


/// \name Member Access Functions 
/// @{

/** @brief returns the `boost::adjacency_list` object for the constructed cone based spanner.
 */
Graph graph();

/** @brief returns the number of cones in this graph. 
 */
const unsigned int& number_of_cones() const;

/** @brief returns the vector of the directions of the rays dividing the plane. 
 *
 *  @return a vector of Direction_2
 */
const std::vector<Direction_2>& directions() const;
	
/// @} 

}; /* end ConeSpanners_2 */

