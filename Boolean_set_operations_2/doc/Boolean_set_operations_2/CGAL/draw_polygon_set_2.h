namespace CGAL {

/*!
\ingroup PkgDrawPolygonSet2

The class `Graphics_scene_options_polygon_set_2` defines data and methods used to tune the way that the cells of a `Polygon_set_2` are considered for drawing or to be added into a graphics scene.

This class is a model of `GraphicsSceneOptionsPolygonSet2`.

\tparam DS a `CGAL::Polygon_set_2`.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.

\cgalModels{GraphicsSceneOptionsPolygonSet2}
*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor>
struct Graphics_scene_options_polygon_set_2: public CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>
{};

/*!
\ingroup PkgDrawPolygonSet2

opens a new window and draws a 2D polygon set. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam PS2 which must be an instantiation of a `CGAL::Polygon_set_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsPolygonSet2` concept.

\param ps2 the polygon set to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class D, class GSOptions>

 void CGAL::draw(const CGAL::Polygon_set_2<T, C, D>& ps2, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class PS2, class GSOptions>
void draw(const PS2& ps2, const GSOptions& gso);

/*!
\ingroup PkgDrawPolygonSet2

A shortcut to `CGAL::draw(ps2, Graphics_scene_options_polygon_set_2{})`.
*/
template<class PS2>
void draw(const PS2& ps2);

/*!
\ingroup PkgDrawPolygonSet2

adds the vertices, edges and faces of `ps2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam PS2 which must be an instantiation of a `CGAL::Polygon_set_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsPolygonSet2` concept.

\param ps2 the polygon set to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class D, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Polygon_set_2<T, C, D>& ps2, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class PS2, class GSOptions>
void add_to_graphics_scene(const PS2& ps2,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawPolygonSet2

A shortcut to `CGAL::add_to_graphics_scene(ps2, gs, Graphics_scene_options_polygon_set_2{})`.
*/
template<class PS2>
void add_to_graphics_scene(const PS2& ps2,
                           CGAL::Graphics_scene& gs);

} /* end namespace CGAL */
