namespace CGAL {

/*!
\ingroup PkgDrawPolyhedron

opens a new window and draws a polyhedron. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam P which must be an instantiation of a `CGAL::Polyhedron_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param p the polyhedron to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class PT, class PI, class HDS, class Alloc, class GSOptions>

 void CGAL::draw(const CGAL::Polyhedron_3<PT, PI, HDS, Alloc>& p, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P, class GSOptions>
void draw(const P& p, const GSOptions& gso);

/*!
\ingroup PkgDrawPolyhedron

A shortcut to `CGAL::draw(p, Graphics_scene_options{})`.
*/
  template<class P>
  void draw(const P& p);

/*!
\ingroup PkgDrawPolyhedron

adds the vertices, edges and faces of `p` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam P which must be an instantiation of a `CGAL::Polyhedron_3<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param p the polyhedron to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class PT, class PI, class HDS, class Alloc, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Polyhedron_3<PT, PI, HDS, Alloc>& p, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P, class GSOptions>
void add_to_graphics_scene(const P& p,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawPolyhedron

A shortcut to `CGAL::add_to_graphics_scene(p, gs, Graphics_scene_options{})`.
*/
template<class P>
void add_to_graphics_scene(const P& p,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */

