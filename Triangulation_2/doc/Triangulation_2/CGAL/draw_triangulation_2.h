namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws a triangulation. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam T2 which must be an instantiation of a `CGAL::Triangulation_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param at2 the triangulation to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class GSOptions>

 void CGAL::draw(const CGAL::Triangulation_2<Gt, Tds>& at2, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
  template<class T2, class GSOptions>
  void draw(const T2& at2, const GSOptions& gso);

/*!
\ingroup PkgDrawTriangulation2

A shortcut to `CGAL::draw(at2, Graphics_scene_options{})`.
*/
  template<class T2>
  void draw(const T2& at2);

/*!
\ingroup PkgDrawTriangulation2

adds the vertices, edges and faces of `at2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam T2 which must be an instantiation of a `CGAL::Triangulation_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param at2 the triangulation to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Triangulation_2<Gt, Tds>& at2, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class T2, class GSOptions>
void add_to_graphics_scene(const T2& at2,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawTriangulation2

A shortcut to `CGAL::add_to_graphics_scene(at2, gs, Graphics_scene_options{})`.
*/
template<class T2>
void add_to_graphics_scene(const T2& at2,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */
