// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2MainClasses

The class `Periodic_2_triangulation_hierarchy_2` implements a
triangulation augmented with a data structure which allows fast point
location queries.

\cgalHeading{Parameters}

It is templated by a parameter which must be instantiated by one of the \cgal periodic triangulation classes. <I>In the current
implementation, only `Periodic_2_Delaunay_triangulation_2` is
supported for `PTr`.</I>

\tparam PTr::Vertex has to be a model of the concept
`Periodic_2TriangulationHierarchyVertexBase_2`.

\tparam PTr::Geom_traits has to be a model of the concept
`Periodic_2DelaunayTriangulationTraits_2`.

\cgalHeading{Inheritance}

`Periodic_2_triangulation_hierarchy_2` offers exactly the same functionalities as `PTr`.
Most of them (point location, insertion, removal \f$ \ldots\f$ ) are overloaded to
improve their efficiency by using the hierarchical structure.

Note that, since the algorithms that are provided are randomized, the
running time of constructing a triangulation with a hierarchy may be
improved when shuffling the data points.

However, the I/O operations are not overloaded. So, writing a
hierarchy into a file will lose the hierarchical structure and reading
it from the file will result in an ordinary triangulation whose
efficiency will be the same as `PTr`.

\cgalHeading{Implementation}

The data structure is a hierarchy of triangulations. The triangulation
at the lowest level is the original triangulation where operations and
point location are to be performed.
Then at each succeeding level, the data structure
stores a triangulation of a small random sample of the vertices
of the triangulation at the preceding level. Point location
is done through a top-down nearest neighbor query.
The nearest neighbor query is first
performed naively in the top level triangulation.
Then, at each following level, the nearest neighbor at that level
is found through a randomized walk performed from
the nearest neighbor found at the preceding level.
Because the number of vertices in each triangulation is only a small
fraction of the number of vertices of the preceding triangulation
the data structure remains small and achieves fast point location
queries on real data.

\sa `CGAL::Periodic_2_triangulation_hierarchy_vertex_base_2`
\sa `CGAL::Periodic_2_Delaunay_triangulation_2`

*/
template< typename PTr >
class Periodic_2_triangulation_hierarchy_2 : public PTr
{
}; /* end Periodic_2_triangulation_hierarchy_2 */
} /* end namespace CGAL */
