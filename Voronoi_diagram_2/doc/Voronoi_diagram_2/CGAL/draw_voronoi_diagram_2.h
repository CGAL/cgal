namespace CGAL {

/*!
\ingroup PkgDrawVoronoiDiagram2

The class `Graphics_scene_options_voronoi_diagram_2` defines data and methods used to tune the way that the cells of a `Voronoi_diagram_2` are considered for drawing or to be added into a graphics scene.

This class is a model of `GraphicsSceneOptionsVoronoiDiagram2`.

\tparam DS a `CGAL::Voronoi_diagram_2`.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.

\cgalModels{GraphicsSceneOptionsVoronoiDiagram2}
*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor>
struct Graphics_scene_options_voronoi_diagram_2: public CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>
{};

/*!
\ingroup PkgDrawVoronoiDiagram2

opens a new window and draws a 2D voronoi diagram. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam VD2 which must be an instantiation of a `CGAL::Voronoi_diagram_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsVoronoiDiagram2` concept.

\param vd2 the voronoi diagram to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class DG, class AT, class AP, class GSOptions>

 void CGAL::draw(const CGAL::Voronoi_diagram_2<DG, AT, AP>& vd2, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class VD2, class GSOptions>
void draw(const VD2& vd2, const GSOptions& gso);

/*!
\ingroup PkgDrawVoronoiDiagram2

A shortcut to `CGAL::draw(vd2, Graphics_scene_options_voronoi_diagram_2{})`.
*/
template<class VD2>
void draw(const VD2& vd2);

/*!
\ingroup PkgDrawVoronoiDiagram2

adds the vertices, edges and faces of `vd2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam VD2 which must be an instantiation of a `CGAL::Voronoi_diagram_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsVoronoiDiagram2` concept.

\param vd2 the voronoi diagram to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class DG, class AT, class AP, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Voronoi_diagram_2<DG, AT, AP>& vd2, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class VD2, class GSOptions>
void add_to_graphics_scene(const VD2& vd2,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawVoronoiDiagram2

A shortcut to `CGAL::add_to_graphics_scene(vd2, gs, Graphics_scene_options_voronoi_diagram_2{})`.
*/
template<class VD2>
void add_to_graphics_scene(const VD2& vd2,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */
