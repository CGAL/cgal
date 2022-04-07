namespace CGAL {

/*!
\ingroup PkgPointSet2NeighborSearch

computes a handle to a vertex `w` of `delau` that is closest to `v`.
If `v` is the only vertex in `delau`, `NULL` is returned.

\cgalHeading{Requirements}

`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from
the Delaunay triangulation data type:
<UL>
<LI>`Dt::Geom_traits`
<LI>`Dt::Point`
<LI>`Dt::Vertex_circulator`
<LI>`Dt::Vertex_handle`
<LI>`Dt::Geom_traits::Compare_distance_2`
</UL>

*/
template<class Dt>
Dt::Vertex_handle nearest_neighbor(const Dt& delau, Dt::Vertex_handle v);

} /* namespace CGAL */

namespace CGAL {

/*!
  \ingroup PkgPointSet2NeighborSearch

computes the `k` nearest neighbors of `p` in `delau`, and places the
handles to the corresponding vertices as a sequence of objects of type
Vertex_handle in a container of value type of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.

The function `nearest_neighbors()` is the function template version of the k nearest
neighbors search on Delaunay triangulations. There are two versions of
this function, one taking a point of the Delaunay triangulation and the
other taking a vertex handle.

\cgalHeading{Requirements}

`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from
the Delaunay triangulation data type:
<UL>
<LI>`Dt::Geom_traits`
<LI>`Dt::Vertex_handle`
<LI>`Dt::Vertex_iterator`
<LI>`Dt::Vertex_circulator`
<LI>`Dt::Vertex`
<LI>`Dt::Face`
<LI>`Dt::Face_handle`
<LI>`Dt::Locate_type`
<LI>`Dt::Point`
<LI>`Dt::Geom_traits::FT`
<LI>`Dt::Geom_traits::Compute_squared_distance_2`
</UL>

*/
template<class Dt, class OutputIterator>
OutputIterator nearest_neighbors(Dt& delau, const Dt::Point& p, Dt::size_type k, OutputIterator res);

/*!
\ingroup PkgPointSet2NeighborSearch

computes the `k` nearest neighbors of `v` (including `v`) in `delau`, and places them as a sequence of objects of type
Vertex_handle in a container of value type of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.

The function `nearest_neighbors()` is the function template version of the k nearest
neighbors search on Delaunay triangulations. There are two versions of
this function, one taking a point of the Delaunay triangulation and the
other taking a vertex handle.

\cgalHeading{Requirements}

`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from
the Delaunay triangulation data type:
<UL>
<LI>`Dt::Geom_traits`
<LI>`Dt::Vertex_handle`
<LI>`Dt::Vertex_iterator`
<LI>`Dt::Vertex_circulator`
<LI>`Dt::Vertex`
<LI>`Dt::Point`
<LI>`Dt::Geom_traits::FT`
<LI>`Dt::Geom_traits::Compute_squared_distance_2`
</UL>

*/
template<class Dt, class OutputIterator>
OutputIterator nearest_neighbors(Dt& delau, Dt::Vertex_handle v, Dt::size_type k, OutputIterator res);

} /* namespace CGAL */
