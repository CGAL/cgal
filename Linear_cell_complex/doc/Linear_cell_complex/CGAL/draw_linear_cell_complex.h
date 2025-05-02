namespace CGAL {

/*!
\ingroup PkgDrawLinearCellComplex

opens a new window and draws a linear cell complex. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam LCC which must be an instantiation of a `CGAL::Linear_cell_complex_for_combinatorial_map<...>` or `CGAL::Linear_cell_complex_for_generalized_map<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param lcc the linear cell complex to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<unsigned int d, unsigned int ad, class T, class I, class M, class R, class S, class GSOptions>

 void CGAL::draw(const CGAL::Linear_cell_complex_base<d, ad, T, I, A, M, R, S>& lcc, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
  template<class LCC, class GSOptions>
  void draw(const LCC& lcc, const GSOptions& gso);

/*!
\ingroup PkgDrawLinearCellComplex

A shortcut to `CGAL::draw(lcc, Graphics_scene_options{})`.
*/
  template<class LCC>
  void draw(const LCC& lcc);

/*!
\ingroup PkgDrawLinearCellComplex

adds the vertices, edges and faces of `lcc` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam LCC which must be an instantiation of a `CGAL::Linear_cell_complex_for_combinatorial_map<...>` or `CGAL::Linear_cell_complex_for_generalized_map<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param lcc the linear cell complex to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<unsigned int d, unsigned int ad, class T, class I, class M, class R, class S, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Linear_cell_complex_base<d, ad, T, I, A, M, R, S>& lcc, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class LCC, class GSOptions>
void add_to_graphics_scene(const LCC& lcc,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawLinearCellComplex

A shortcut to `CGAL::add_to_graphics_scene(lcc, gs, Graphics_scene_options{})`.
*/
template<class LCC>
void add_to_graphics_scene(const LCC& lcc,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */

