namespace CGAL {

/*!
\defgroup PkgInterpolationNaturalNeighborCoordinates2 CGAL::natural_neighbor_coordinates_2()
\ingroup PkgInterpolation2NatNeighbor

The functions `natural_neighbor_coordinates_2()` compute natural neighbor coordinates, also
called Sibson's coordinates, for `2D` points provided a two-dimensional triangulation
and a query point in the convex hull of the vertices of the triangulation.

\cgalHeading{Implementation}

This function computes the area of the sub-cells stolen from the Voronoi cells
of the points in `dt` when inserting `p`. The total area of the Voronoi cell of `p`
is also computed and returned by the function, as second value of a triple.

\attention If `p` lies outside the convex hull, the coordinate values cannot be computed,
and the third value of the resulting triple is set to `false` (the other values
of the triple are then meaningless).

\cgalHeading{Output Format}

The return type is identical for all overloads of `CGAL::natural_neighbor_coordinates_2()`:
it is `CGAL::Triple<CoordinateOutputIterator, Dt::Geom_traits::FT, bool >`

Natural neighbor coordinates are output in the first value of the triple,
using an output iterator (see the concept `OutputIterator`).
Internally, these coordinates are associated to some of the vertices of the triangulation,
and a natural value type for the output iterator would thus be
`std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`.
However, to allow flexibility in the format of the output, a functor passed as argument
can be used to format the output as desired: this functor must have argument type
`std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`
and its result type and the value type of the output iterator must be identical.

For example, we can use a functor of type
\link CGAL::Identity
CGAL::Identity<std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT> > \endlink
and an output iterator with value type `std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`
or we could extract the bare point and use a value type: `std::pair<Dt::Point, Dt::Geom_traits::FT>`,
or simply the coordinates, etc.

This functor can be omitted and will default, for backward compatibility reasons,
to a functor that extracts the point from the vertex, and the output iterator
must then have value type `std::pair<Dt::Point, Dt::Geom_traits::FT>`.

\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa `PkgInterpolationSurfaceNeighborCoordinates3`
\sa `PkgInterpolation2Interpolation`
*/
/// @{

/*!
Computes the natural neighbor coordinates for `p` with respect to the points
in the two-dimensional Delaunay triangulation `dt`.

\tparam Dt must be of type `Delaunay_triangulation_2<Traits, Tds>`.
        `Traits` must be a model of the concepts `DelaunayTriangulationTraits_2` and `PolygonTraits_2`.
\tparam CoordinateOutputIterator must be a model of `OutputIterator` and have the value type `OutputFunctor::result_type`.
        The output computed by the function is placed starting at `out`.
\tparam OutputFunctor must be a functor with argument type `std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`.
        See \link PkgInterpolationNaturalNeighborCoordinates2 above \endlink for a detailed explanation on the usage of `OutputFunctor`.

\param dt is the Delaunay triangulation.
\param p is the query point.
\param out is an object of type `CoordinateOutputIterator`.
\param fct is an object of type `OutputFunctor`.
\param start is an optional argument that is used as a hint of where the locate process has to start its search.

\return A triple consisting of:
- a sequence of objects of types `OutputFunctor::result_type`, starting at `out`.
- the normalization factor of the coordinates.
- a Boolean value which is set to `true` if the coordinate computation was successful,
and `false` otherwise.
*/
template < class Dt, class CoordinateOutputIterator, class OutputFunctor >
CGAL::Triple<CoordinateOutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               CoordinateOutputIterator out,
                               OutputFunctor fct,
                               typename Dt::Face_handle start = typename Dt::Face_handle());

/*!
Computes the natural neighbor coordinates for `p` with respect to the points
in the two-dimensional Delaunay triangulation `dt`.

The iterator range `[hole_begin, hole_end)` determines the boundary edges
of the conflict zone of `p` in the triangulation. It can, for example, be computed
using the function: \link Delaunay_triangulation_2::get_boundary_of_conflicts()
`dt.get_boundary_of_conflicts(p, std::back_inserter(hole), start)`\endlink.

\cgalHeading{Requirements}
Same as above.
*/
template < class Dt, class CoordinateOutputIterator, class OutputFunctor, class EdgeIterator >
CGAL::Triple< CoordinateOutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               CoordinateOutputIterator out,
                               OutputFunctor fct,
                               EdgeIterator hole_begin, EdgeIterator hole_end);

/*!
Computes the natural neighbor coordinates of the point `vh->point()`
with respect to the vertices of `dt` excluding `vh->point()`.

\cgalHeading{Requirements}
Same as above.
*/
template < class Dt, class CoordinateOutputIterator, class OutputFunctor >
CGAL::Triple< CoordinateOutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               typename Dt::Vertex_handle vh,
                               CoordinateOutputIterator out,
                               OutputFunctor fct);

/// @}

} /* namespace CGAL */
