
/*!
\ingroup PkgTriangulation2Concepts
\cgalconcept

The concept `ConstrainedDelaunayTriangulationTraits_2` defines the requirements for the geometric 
traits class of a constrained Delaunay triangulation 
that supports intersections of input constraints. 
This is the case 
when the template parameter `Itag` 
of ( `Constrained_Delaunay_triangulation_2<Traits,Tds,Itag>`) 
is instantiated 
by one of the tag classes `Exact_intersections_tag` or 
`Exact_predicates_tag`). 
The concept `ConstrainedDelaunayTriangulationTraits_2` refines both the concept 
`DelaunayTriangulationTraits_2` and the concept 
`ConstrainedTriangulationTraits_2`. 

\refines ::DelaunayTriangulationTraits_2 
\refines ::ConstrainedTriangulationTraits_2 

\hasModel All \cgal Kernels
\hasModel `CGAL::Triangulation_euclidean_traits_xy<K>` 
\hasModel `CGAL::Triangulation_euclidean_traits_yz<K>` 
\hasModel `CGAL::Triangulation_euclidean_traits_xz<K>` 

\sa `TriangulationTraits_2` 
\sa `ConstrainedTriangulationTraits_2` 
\sa `CGAL::Constrained_triangulation_2<Gt,Tds,Itag>` 
\sa `CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>` 

*/

class ConstrainedDelaunayTriangulationTraits_2 {
public:

/// @}

}; /* end ConstrainedDelaunayTriangulationTraits_2 */

