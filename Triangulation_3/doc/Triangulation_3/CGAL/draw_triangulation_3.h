namespace CGAL {

/*!
\ingroup PkgDrawTriangulation3

opens a new window and draws a 3D triangulation. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam T3 which must be an instantiation of a `CGAL::Triangulation_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param at3 the triangulation to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class LDS, class GSOptions>

 void CGAL::draw(const CGAL::Triangulation_3<Gt, Tds, LDS>& at3, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class T3, class GSOptions>
void draw(const T3& at3, const GSOptions& gso);

/*!
\ingroup PkgDrawTriangulation3
 A shortcut to `CGAL::draw(at3, Graphics_scene_options{})`.
*/
  template<class T3>
  void draw(const T3& at3);

/*!
\ingroup PkgDrawTriangulation3

adds the vertices, edges and faces of `at3` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam T3 which must be an instantiation of a `CGAL::Triangulation_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param at3 the triangulation to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class LDS, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Triangulation_3<Gt, Tds, LDS>&  at3, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class T3, class GSOptions>
void add_to_graphics_scene(const T3& at3,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawTriangulation3

A shortcut to `CGAL::add_to_graphics_scene(at3, gs, Graphics_scene_options{})`.
*/
template<class T3>
void add_to_graphics_scene(const T3& at3,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */

