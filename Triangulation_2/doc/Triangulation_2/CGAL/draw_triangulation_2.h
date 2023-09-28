namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws a triangulation. Parameters of the drawing are taken from the optional cell parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with
`CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam T2 a triangulation class derived from `Triangulation_2`
\param at2 the triangulation to draw.
\param cp the cell parameters.

*/
template<class T2>
void draw(const T2& at2, const CGAL::Graphics_scene_options<T2,
                  typename T2::Vertex_handle,
                  typename T2::Finite_edges_iterator,
                  typename T2::Finite_faces_iterator>& cp=default);

/*!
\ingroup PkgDrawTriangulation2

adds the vertices, edges and faces of `at2` into the given graphic storage `gs`. Parameters of the cells are taken from the optional cell parameters `cp` . Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam T2 a triangulation class derived from `Triangulation_2`
\param at2 the triangulation to draw.
\param gs the graphic storage to fill.
\param cp the cell parameters.

*/
template<class T2>
void add_in_graphic_storage(const T2& at2,
                            CGAL::Graphics_scene<BufferType>& gs,
                            const CGAL::Graphics_scene_options<T2,
                            typename T2::Vertex_handle,
                            typename T2::Finite_edges_iterator,
                            typename T2::Finite_faces_iterator>& cp=default);

} /* namespace CGAL */
