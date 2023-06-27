namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws a triangulation. Parameters of the drawing are taken from the optional drawing functor parameter.

A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with
`CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam T2 a triangulation class derived from `Triangulation_2`
\param at2 the triangulation to draw.
\param df the drawing functor.

*/
template<class T2>
void draw(const T2& at2, const CGAL::Drawing_functor<T2,
                  typename T2::Vertex_handle,
                  typename T2::Finite_edges_iterator,
                  typename T2::Finite_faces_iterator>& df=default);

/*!
\ingroup PkgDrawTriangulation2

adds the vertices, edges and faces of at2 into the given graphic storage gs. Parameters of the cells are taken from the optional drawing functor parameter. Note that gs is not clear before to be filled (to enable to draw several data structures in a same basic viewer).

\tparam T2 a triangulation class derived from `Triangulation_2`
\param at2 the triangulation to draw.
\param gs the graphic storage to fill.
\param df the drawing functor.

*/
template<class T2>
void add_in_graphic_storage(const T2& at2,
                            CGAL::Graphic_storage<BufferType>& gs,
                            const CGAL::Drawing_functor<T2,
                            typename T2::Vertex_handle,
                            typename T2::Finite_edges_iterator,
                            typename T2::Finite_faces_iterator>& df=default);

} /* namespace CGAL */
