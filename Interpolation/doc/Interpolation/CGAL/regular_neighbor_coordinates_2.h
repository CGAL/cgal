namespace CGAL {

/*!
\defgroup PkgInterpolationRegularNeighborCoordinates2 CGAL::regular_neighbor_coordinates_2()
\ingroup PkgInterpolation2NatNeighbor

The functions `regular_neighbor_coordinates_2()` compute natural neighbor coordinates, also 
called Sibson's coordinates, for weighted `2D` points provided a 
two-dimensional regular triangulation and a (weighted) query point 
inside the convex hull of the vertices of the triangulation. We call these 
coordinates regular neighbor coordinates. 

\cgalHeading{Requirements}

<OL> 
<LI>`Rt` are equivalent to the class 
`Regular_triangulation_2<Traits, Tds>`. 
<LI>The traits class `Traits` of `Rt` is a model of the 
concept `RegularTriangulationTraits_2`. It provides the number 
type `FT` which is a model for `FieldNumberType` and it must 
meet the requirements for the traits class of the 
`polygon_area_2()` function. A model of this traits class is 
`Regular_triangulation_euclidean_traits_2<K, Weight>`. 
<LI>The value type of `OutputIterator` is equivalent to 
`std::pair<Rt::Weighted_point, Rt::Geom_traits::FT>`, i.e.\ a pair
associating a point and its regular neighbor coordinate. 
</OL> 

\cgalHeading{Implementation}

This function computes the areas stolen from the 
Voronoi cells of points in `rt` by the insertion of `p`. The 
total area of the Voronoi cell of `p` is also computed and 
returned by the function. If `p` lies outside the convex hull, the 
coordinate values cannot be computed and the third value of the result 
triple is set to `false`. 

\sa `PkgInterpolationNaturalNeighborCoordinates2`

*/
/// @{

/*!
computes the regular neighbor coordinates for `p` with respect
to the weighted points in the two-dimensional regular triangulation
`rt`. 

\tparam Rt must be a `Regular_triangulation_2<Traits, Tds>`.
\tparam OutputIterator must have the value type
`std::pair<Rt::Weighted_point, Rt::Geom_traits::FT>`. The sequence of
point/coordinate pairs that is computed by the function is placed
starting at `out`. 

The function returns a triple with an
iterator that is placed past-the-end of the resulting sequence of
point/coordinate pairs, the normalization factor of the coordinates
and a Boolean value which is set to `true`, iff the coordinate
computation was successful, i.e., if `p` lies inside the
convex hull of the points in `rt`. 
*/
template < class Rt, class OutputIterator > CGAL::Triple<
OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt, const typename
Rt::Weighted_point& p, OutputIterator out, typename
Rt::Face_handle start = typename Rt::Face_handle());

/*!
The same as above. The iterator range `[hole_begin, hole_end)` determines 
the boundary edges of the conflict zone of `p` in the triangulation `rt`.  
\link  Regular_triangulation_2::hidden_vertices_begin() `rt.hidden_vertices_begin()`\endlink and 
\link  Regular_triangulation_2::hidden_vertices_end() `rt.hidden_vertices_end()`\endlink
determines the iterator range over the hidden vertices of the conflict
zone of `p` in`rt`. It is the result of the function
\link Regular_triangulation_2::get_boundary_of_conflicts() `rt.get_boundary_of_conflicts(p,std::back_inserter(hole), std::back_inserter(hidden_vertices), start)`\endlink.
*/
template <class Rt, class OutputIterator,
class EdgeIterator, class VertexIterator > CGAL::Triple<
OutputIterator, typename Traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt, const typename
Traits::Weighted_point& p, OutputIterator out, EdgeIterator
hole_begin, EdgeIterator hole_end, VertexIterator
hidden_vertices_begin, VertexIterator hidden_vertices_end);

/*!
computes the regular neighbor coordinates of the point
`vh->point()` with respect to the vertices of `rt` excluding
`vh->point()`. The same as above for the remaining parameters.
*/
template <class Rt, class OutputIterator>
CGAL::Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt, typename
Rt::Vertex_handle vh, OutputIterator out);

/// @}

} /* namespace CGAL */

