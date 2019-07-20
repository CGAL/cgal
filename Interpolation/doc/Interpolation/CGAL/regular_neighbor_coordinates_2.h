namespace CGAL {

/*!
\defgroup PkgInterpolationRegularNeighborCoordinates2 CGAL::regular_neighbor_coordinates_2()
\ingroup PkgInterpolation2NatNeighbor

The functions `regular_neighbor_coordinates_2()` compute natural neighbor coordinates, also
called Sibson's coordinates, for weighted `2D` points provided a two-dimensional regular triangulation
and a (weighted) query point inside the convex hull of the vertices of the triangulation.
We call these coordinates <em>regular neighbor coordinates</em>.

\cgalHeading{Implementation}

This function computes the areas stolen from the power cells of points
in the regular triangulation `rt` by the insertion of the query point `p`.
The total area of the power cell of `p` is also computed and returned by the function.

\attention If `p` lies outside the convex hull, the coordinate values cannot be
computed and the third value of the result triple is set to `false`.

\cgalHeading{Output Format}

The return type is identical for all overloads of `CGAL::regular_neighbor_coordinates_2()`:
it is `CGAL::Triple<CoordinateOutputIterator, Rt::Geom_traits::FT, bool >`

Regular neighbor coordinates are output in the first value of the triple,
using an output iterator (see the concept `OutputIterator`).
Internally, these coordinates are associated to some of the vertices of the triangulation,
and a natural value type for the output iterator would thus be
`std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`.
However, to allow flexibility in the format of the output, a functor passed as argument
can be used to format the output as desired: this functor must have argument type
`std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
and its result type and the value type of the output iterator must be identical.

For example, we can use a functor of type
\link CGAL::Identity
CGAL::Identity<std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT> > \endlink
and an output iterator with value type `std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
or we could extract the bare point and use a value type: `std::pair<Rt::Weighted_point, Rt::Geom_traits::FT>`,
or simply the coordinates, etc.

This functor can be omitted and will default, for backward compatibility reasons,
to a functor that extracts the point from the vertex, and the output iterator
must then have value type `std::pair<Rt::Weighted_point, Rt::Geom_traits::FT>`.

\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationSurfaceNeighborCoordinates3`
\sa `PkgInterpolation2Interpolation`
*/
/// @{

/*!
Computes the regular neighbor coordinates for `p` with respect
to the weighted points in the two-dimensional regular triangulation `rt`.

\tparam Rt must be a `Regular_triangulation_2<Traits, Tds>`.
        `Traits` must be a model of the concepts `RegularTriangulationTraits_2` and `PolygonTraits_2`.
\tparam CoordinateOutputIterator must be a model of `OutputIterator` and have the value type `OutputFunctor::result_type`.
        The output computed by the function is placed starting at `out`.
\tparam OutputFunctor must be a functor with argument type `std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`.

See \link PkgInterpolationRegularNeighborCoordinates2 above \endlink
for a detailed explanation on the usage of `OutputFunctor`.

\param rt is the regular triangulation.
\param p is the query point.
\param out is an object of type `CoordinateOutputIterator`.
\param fct is an object of type `OutputFunctor`.
\param start is an optional argument that is used as a hint of where the locate process has to start its search.

\return A triple consisting of:
- a sequence of objects of types `OutputFunctor::result_type`, starting at `out`.
- the normalization factor of the coordinates.
- a Boolean value which is set to `true` if the coordinate computation was successful,
and `false` otherwise.

\warning If the weight of `p` is so small that the point `p` would not have
any power cell if it were inserted in the power diagram, then the resulting triple
will be `(out, 0, true)` with no additional entry in `out` (compared to its state
in input).
*/
template < class Rt, class CoordinateOutputIterator, class OutputFunctor >
CGAL::Triple<CoordinateOutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               CoordinateOutputIterator out,
                               OutputFunctor fct,
                               typename Rt::Face_handle start = typename Rt::Face_handle());

/*!
Computes the regular neighbor coordinates for `p` with respect
to the weighted points in the two-dimensional regular triangulation `rt`.

The iterator range `[hole_begin, hole_end)` determines the boundary edges
of the conflict zone of `p` in the triangulation `rt`.
\link Regular_triangulation_2::hidden_vertices_begin() `rt.hidden_vertices_begin()`\endlink and
\link Regular_triangulation_2::hidden_vertices_end() `rt.hidden_vertices_end()`\endlink
determines the iterator range over the hidden vertices of the conflict
zone of `p` in`rt`. Those ranges can, for example, be computed using the function
\link Regular_triangulation_2::get_boundary_of_conflicts_and_hidden_vertices()
`rt.get_boundary_of_conflicts_and_hidden_vertices(p,std::back_inserter(hole),
std::back_inserter(hidden_vertices), start)`\endlink.

\cgalHeading{Requirements}
Same as above.
*/
template <class Rt, class CoordinateOutputIterator, class OutputFunctor, class EdgeIterator, class VertexIterator >
CGAL::Triple< CoordinateOutputIterator, typename Traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Traits::Weighted_point& p,
                               CoordinateOutputIterator out,
                               OutputFunctor fct,
                               EdgeIterator hole_begin, EdgeIterator hole_end,
                               VertexIterator hidden_vertices_begin, VertexIterator hidden_vertices_end);

/*!
Computes the regular neighbor coordinates of the point `vh->point()` with respect
to the vertices of `rt` excluding `vh->point()`.

\cgalHeading{Requirements}
Same as above.
*/
template < class Rt, class CoordinateOutputIterator, class OutputFunctor >
CGAL::Triple< CoordinateOutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               typename Rt::Vertex_handle vh,
                               CoordinateOutputIterator out,
                               OutputFunctor fct);

/// @}

} /* namespace CGAL */
