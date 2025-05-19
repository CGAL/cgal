namespace CGAL {

/*!
\ingroup PkgDrawNef3

opens a new window and draws a nef polyhedron. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam NP3 which must be an instantiation of a `CGAL::Nef_polyhedron_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param np3 the nef polyhedron to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class K, class I, class M, class GSOptions>

 void CGAL::draw(const CGAL::Nef_polyhedron_3<K, I, M>& np3, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
  template<class NP3, class GSOptions>
  void draw(const NP3& np3, const GSOptions& gso);

/*!
\ingroup PkgDrawNef3

A shortcut to `CGAL::draw(np3, Graphics_scene_options{})`.
*/
  template<class NP3>
  void draw(const NP3& np3);

/*!
\ingroup PkgDrawNef3

adds the vertices, edges and faces of `np3` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam NP3 which must be an instantiation of a `CGAL::Nef_polyhedron_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param np3 the nef polyhedron to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class K, class I, class M, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Nef_polyhedron_3<K, I, M>& np3, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class NP3, class GSOptions>
void add_to_graphics_scene(const NP3& np3,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawNef3

A shortcut to `CGAL::add_to_graphics_scene(np3, gs, Graphics_scene_options{})`.
*/
template<class NP3>
void add_to_graphics_scene(const NP3& np3,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */
