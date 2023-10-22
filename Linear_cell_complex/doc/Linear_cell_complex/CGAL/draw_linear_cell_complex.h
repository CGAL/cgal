namespace CGAL {

/*!
\ingroup PkgDrawLinearCellComplex

opens a new window and draws `alcc`, a model of the `LinearCellComplex` concept. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam LCC a model of the `LinearCellComplex` concept.
\tparam GSOptions a class having the same methods than `Graphics_scene_options`: `Graphics_scene_options` by default.

\param alcc the linear cell complex to draw.
\param gs_options graphics scene options.

*/
  template<class LCC, class GSOptions>
  void draw(const LCC& alcc, const GSOptions& gs_options=GSOptions());

/*!
\ingroup PkgDrawLinearCellComplex

Add in the given graphics scene the elements of alcc.

\tparam LCC a model of the `LinearCellComplex` concept.
\tparam BufferType the number type used for point coordinates: `float` by default.
\tparam GSOptions a class having the same methods than `Graphics_scene_options`: `Graphics_scene_options` by default.

\param alcc the linear cell complex to draw.
\param graphics_scene the graphics scene.
\param gs_options graphics scene options.

*/
  template<class LCC, class BufferType, class GSOptions>
  void add_in_graphics_scene(const LCC& alcc, CGAL::Graphics_scene& graphics_scene,
                             const GSOptions& gs_options=GSOptions());

} /* namespace CGAL */

