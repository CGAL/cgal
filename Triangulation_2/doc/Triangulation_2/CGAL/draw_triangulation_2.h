namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws a triangulation. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with
`CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam T2 a triangulation class derived from `Triangulation_2`
\param at2 the triangulation to draw.
\param gso the graphics scene options parameter, `Graphics_scene_options()` by default.

*/
template<class T2>
void draw(const T2& at2,
          const CGAL::Graphics_scene_options<T2,
          typename T2::Vertex_handle,
          typename T2::Finite_edges_iterator,
          typename T2::Finite_faces_iterator>& gso=default);

/*!
\ingroup PkgDrawTriangulation2

adds the vertices, edges and faces of `at2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso` . Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam T2 a triangulation class derived from `Triangulation_2`
\tparam BufferType the number type used for point coordinates: `float` by default.

\param at2 the triangulation to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter, `Graphics_scene_options()` by default.

*/
template<class T2, typename BufferType=float>
void add_to_graphics_scene(const T2& at2,
                           CGAL::Graphics_scene<BufferType>& gs,
                           const CGAL::Graphics_scene_options<T2,
                           typename T2::Vertex_handle,
                           typename T2::Finite_edges_iterator,
                           typename T2::Finite_faces_iterator>& gso=default);

} /* namespace CGAL */
