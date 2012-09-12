namespace CGAL {

/*!
\defgroup natural_neighbor_coordinates_2 natural_neighbor_coordinates_2
\ingroup PkgInterpolation2NatNeighbor

The function `natural_neighbor_coordinates_2` computes natural neighbor coordinates, also 
called Sibson's coordinates, for \f$ 2D\f$ points provided a two-dimensional 
triangulation and a query point in the convex hull of the vertices 
of the triangulation. 

Requirements 
-------------- 

<OL> 
<LI>`Dt` are equivalent to the class 
`Delaunay_triangulation_2<Traits, Tds>`. 
<LI>The traits class `Traits` of `Dt` is a model of the 
concept `DelaunayTriangulationTraits_2`. 
Only the following members of this traits class are used: 
<UL> 
<LI>`Construct_circumcenter_2` 
<LI>`FT` 
<LI>`Point_2` 
<LI>`construct_circumcenter_2_object` 
<DT><B></B><DD> Additionally, `Traits` must meet the requirements for 
the traits class of the `polygon_area_2` function. 
</UL> 
<LI>`OutputIterator::value_type` is equivalent to 
`std::pair<Dt::Point_2, Dt::Geom_traits::FT>`, i.e. a pair 
associating a point and its natural neighbor coordinate. 
</OL> 

Implementation 
-------------- 

This function computes the area of the sub-cells 
stolen from the Voronoi cells of the points in `dt` when inserting 
`p`. The total area of the Voronoi cell of `p` is also 
computed and returned by the function. If `p` lies outside the 
convex hull, the coordinate values cannot be computed and the third 
value of the result triple is set to `false`. 


\sa CGAL::linear_interpolation 
\sa CGAL::sibson_c1_interpolation 
\sa CGAL::surface_neighbor_coordinates_3 
\sa CGAL::regular_neighbor_coordinates_2 

*/
/// @{

/*!
computes the natural neighbor coordinates for `p` with respect to the
points in the two-dimensional Delaunay triangulation `dt`. The
template class `Dt` should be of type
`Delaunay_triangulation_2<Traits, Tds>`. The value type of the
`OutputIterator` is a pair of `Dt::Point_2` and the coordinate value
of type `Dt::Geom_traits::FT`. The sequence of point/coordinate pairs
that is computed by the function is placed starting at `out`. The
function returns a triple with an iterator that is placed past-the-end
of the resulting sequence of point/coordinate pairs, the normalization
factor of the coordinates and a Boolean value which is set to true iff
the coordinate computation was successful.
*/
template < class Dt, class OutputIterator > CGAL::Triple<
OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(
  const Dt& dt, const typename Dt::Geom_traits::Point_2& p, 
  OutputIterator out, typename Dt::Face_handle start = typename Dt::Face_handle());

/*!
The same as above. `hole_begin` and `hole_end` determines the
iterator range over the boundary edges of the conflict zone of `p` in
the triangulation. It is the result of the function
`T.get_boundary_of_conflicts(p,std::back_inserter(hole), start)`, see
`Delaunay_triangulation_2<Traits, Tds>`.
*/
template <class Dt, class OutputIterator,
class EdgeIterator > CGAL::Triple< OutputIterator, typename Dt::Geom_traits::FT, 
bool > natural_neighbor_coordinates_2(
  const Dt& dt, const typename Dt::Geom_traits::Point_2& p, 
  OutputIterator out, EdgeIterator hole_begin, EdgeIterator hole_end);

/*!
This function computes the natural neighbor coordinates of the point
`vh->point()` with respect to the vertices of `dt` excluding
`vh->point()`. The same as above for the remaining parameters.
*/
template <class Dt, class OutputIterator> 
CGAL::Triple< OutputIterator, typename Dt::Geom_traits::FT, bool > 
natural_neighbor_coordinates_2(const Dt& dt, typename Dt::Vertex_handle vh, OutputIterator out);

/// @}

} /* namespace CGAL */

