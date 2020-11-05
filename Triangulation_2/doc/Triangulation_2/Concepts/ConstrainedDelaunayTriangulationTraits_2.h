
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The concept `ConstrainedDelaunayTriangulationTraits_2` defines the requirements for the geometric
traits class of a constrained Delaunay triangulation
that supports intersections of input constraints.
This is the case
when the template parameter `Itag`
of `CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds,Itag>` is instantiated
by one of the tag classes `CGAL::Exact_intersections_tag` or
`CGAL::Exact_predicates_tag`.
The concept `ConstrainedDelaunayTriangulationTraits_2` refines both the concept
`DelaunayTriangulationTraits_2` and the concept
`ConstrainedTriangulationTraits_2`.

\cgalRefines `DelaunayTriangulationTraits_2`
\cgalRefines `ConstrainedTriangulationTraits_2`

\cgalHasModel All \cgal Kernels
\cgalHasModel `CGAL::Projection_traits_xy_3<K>`
\cgalHasModel `CGAL::Projection_traits_yz_3<K>`
\cgalHasModel `CGAL::Projection_traits_xz_3<K>`


\sa `TriangulationTraits_2`
\sa `ConstrainedTriangulationTraits_2`
\sa `CGAL::Constrained_triangulation_2<Gt,Tds,Itag>`
\sa `CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>`

*/

class ConstrainedDelaunayTriangulationTraits_2 {
public:


}; /* end ConstrainedDelaunayTriangulationTraits_2 */

