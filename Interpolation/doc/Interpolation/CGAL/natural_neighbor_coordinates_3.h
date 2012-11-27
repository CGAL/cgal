namespace CGAL {

/*!
\defgroup natural_neighbor_coordinates_3 natural_neighbor_coordinates_3
\ingroup PkgInterpolation2NatNeighbor

Given a 3D point `p` and a 3D Delaunay triangulation `dt`, 
 this function calculates the natural neighbors and coordinates of `p` with regard of `dt`.
 
\tparam OutputIterator must have value type
       `std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`

Result : 
- An `OutputIterator` providing natural neighbors `p_i` of `p` with unnormalized coordinates `a_i` associated to them
- The normalizing coefficient (sum over i of the a_i)
- A boolean specifying whether the calculation has succeeded or not
*/

/// @{

/*!
 */
template <class Dt, class OutputIterator>
Triple< OutputIterator, 
	typename Dt::Geom_traits::FT,
	bool >
laplace_natural_neighbor_coordinates_3(const Dt& dt,
				       const typename Dt::Geom_traits::Point_3& p,
				       OutputIterator nn_out, typename Dt::Geom_traits::FT & norm_coeff,
				       const typename Dt::Cell_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Cell_handle());

  /*!
   */
template <class Dt, class OutputIterator>
Triple< OutputIterator,
	typename Dt::Geom_traits::FT,
	bool >
sibson_natural_neighbor_coordinates_3(const Dt& dt,
				      const typename Dt::Geom_traits::Point_3& p,
				      OutputIterator nn_out, typename Dt::Geom_traits::FT & norm_coeff,
				      const typename Dt::Cell_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Cell_handle());

  /// @}
}
