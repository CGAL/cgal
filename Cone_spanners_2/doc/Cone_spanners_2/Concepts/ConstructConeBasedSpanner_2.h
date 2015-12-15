/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The concept `ConstructConeBasedSpanner_2` describes the set of requirements to be fulfilled by any functor class that 
constructs a cone based spanner, such as a Theta graph and a Yao graph, etc. 

\cgalHasModel `CGAL::Construct_theta_graph_2`
\cgalHasModel `CGAL::Construct_yao_graph_2`

*/
template <typename ConeBasedSpannerTraits_2, typename Graph_>
class ConstructConeBasedSpanner_2 {
public:

/// \name Types 
/// @{

/*! The CGAL kernel type used by the functor. Its requirements are described in
    the concept `ConeBasedSpannerTraits_2`.
*/
typedef ConeBasedSpannerTraits_2      Kernel_type;

/*! The graph type to store the constructed cone based spanner. 
    This package requires the use of the `boost::adjacency_list` in the Boost Graph Library as 
	the graph type. Note that there are seven template parameters for
    `boost::adjacency_list`: `OutEdgeList`, `VertexList`, `Directed`, `VertexProperties`, `EdgeProperties`,
    `GraphProperties`, `EdgeList`, of which we require the `VertexProperties` be `CGAL::Point_2`,
    while other parameters can be chosen freely. Here the CGAL kernel type used by `CGAL::Point_2` is
	determined by the first template parameter `ConeBasedSpannerTraits_2`, and we pass `CGAL::Point_2` directly 
	to `adjacency_list` as bundled properties because this makes our implementation much more 
	straightforward than using property maps.
	For detailed information about bundled properties, please refer to
	http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/bundles.html.
	If more properties for vertices are needed, they can be added later as external properties using 
	property maps.
*/
typedef Graph_                           Graph_type;
/*! The directon type.
*/
  typedef ConeBasedSpannerTraits_2::Direction_2      Direction_2;

/// @}


/// \name Creation 
/// @{

/*! \brief Constructor.

   \param k     Number of cones to divide the plane
   \param initial_direction  A direction denoting one of the rays dividing the
                  cones. This allows arbitary rotations of the rays that divide
	              the plane.  (default: positive x-axis)
*/
Construct_spanner_2(unsigned int k, Direction_2 initial_direction = Direction_2(1,0) );

/// @} 

/// \name Operator 
/// @{

/*! \brief Function  operator to construct a cone based spanner.
 *
 * \tparam PointInputIterator This template parameter is to give the application developer freedom 
 *          in choosing whichever iterator he will use to input the coordinates of the points. 
 *          If omitted, the type of the Iterator can be inferred from the arguments passed to the 
 *          operator().
 *
 * \param[in] start An iterator pointing to the first point (vertex).
 * \param[in] end   Past-the-end iterator.
 * \param[out] g   The constructed graph object.
 */
template <typename PointInputIterator>
Graph_& operator()(const PointInputIterator& start,
				   const PointInputIterator& end,
				   Graph_& g);

	
/// @} 

/// \name Member Access Functions 
/// @{

/*! \brief returns the number of cones in this graph. 
 */
const unsigned int number_of_cones() const;

/*! \brief returns the vector of the directions of the rays dividing the plane. 
 */
const std::vector<Direction_2>& directions() const;
	
/// @} 

}; /* end ConstructConeBasedSpanner_2 */

