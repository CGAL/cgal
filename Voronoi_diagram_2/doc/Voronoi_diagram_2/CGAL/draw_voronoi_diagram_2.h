namespace CGAL {

/*!
\ingroup PkgDrawVoronoiDiagram2

opens a new window and draws `av2`, the `Voronoi_diagram_2` constructed from a Delaunay Graph which is a model of `DelaunayGraph_2` concept.
The class `Voronoi_diagram_2` provides an adaptor to view a triangulated Delaunay graph as their dual subdivision, the
Voronoi diagram. A call to this function is blocking, that is the program continues as soon as the user closes the window.
This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam V2 a model of the `AdaptationTraits_2` concept.
\param av2 the voronoi diagram to draw.

*/
template<class V2>
void draw(const V2& av2);

} /* namespace CGAL */
