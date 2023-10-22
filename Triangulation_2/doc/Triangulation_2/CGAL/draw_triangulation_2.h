namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws a triangulation. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with
`CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam Gt
\tparam Tds
\tparam GSOptions a model of `GraphicsSceneOptions` concept, `Graphics_scene_options` by defafult.

\param at2 the triangulation to draw.
\param gso the graphics scene options parameter, `GSOptions()` by default.

*/
  template<class Gt, class Tds, class GSOptions>
  void draw(const CGAL::Triangulation_2<Gt, Tds>& at2,
            const GSOptions& gso={});

/*!
\ingroup PkgDrawTriangulation2

adds the vertices, edges and faces of `at2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso` . Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam Gt
\tparam Tds
\tparam GSOptions a model of `GraphicsSceneOptions` concept, `Graphics_scene_options` by defafult.

\param at2 the triangulation to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter, `GSOptions()` by default.
*/
template<class Gt, class Tds, class GSOptions>
void add_to_graphics_scene(const CGAL::Triangulation_2<Gt, Tds>& at2,
                           CGAL::Graphics_scene<BufferType>& gs,
                           const GSOptions& gso={});

} /* namespace CGAL */
