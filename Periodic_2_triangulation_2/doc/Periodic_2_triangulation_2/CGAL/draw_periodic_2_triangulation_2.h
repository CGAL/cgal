namespace CGAL {

/*!
\ingroup PkgDrawPeriodic2Triangulation2

The class `Graphics_scene_options_periodic_2_triangulation_2` defines data and methods used to tune the way that the cells of a `Periodic_2_triangulation_2` are considered for drawing or to be added into a graphics scene.

This class is a model of `GraphicsSceneOptionsPeriodic2Triangulation2`.

\tparam DS a `CGAL::Periodic_2_triangulation_2`.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.

\cgalModels{GraphicsSceneOptionsPeriodic2Triangulation2}
*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor>
struct Graphics_scene_options_periodic_2_triangulation_2: public CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>
{};

/*!
\ingroup PkgDrawPeriodic2Triangulation2

opens a new window and draws a periodic 2D triangulation. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam P2T2 which must be an instantiation of a `CGAL::Periodic_2_triangulation_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsPeriodic2Triangulation2` concept.

\param p2t2 the periodic triangulation to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class GSOptions>

 void CGAL::draw(const CGAL::Periodic_2_triangulation_2<Gt, Tds>& p2t2, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P2T2, class GSOptions>
void draw(const P2T2& p2t2, const GSOptions& gso);

/*!
\ingroup PkgDrawPeriodic2Triangulation2

A shortcut to `CGAL::draw(p2t2, Graphics_scene_options_periodic_2_triangulation_2{})`.
*/
template<class P2T2>
void draw(const P2T2& p2t2);

/*!
\ingroup PkgDrawPeriodic2Triangulation2

adds the vertices, edges and faces of `p2t2` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam P2T2 which must be an instantiation of a `CGAL::Periodic_2_triangulation_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptionsPeriodic2Triangulation2` concept.

\param p2t2 the periodic triangulation to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Periodic_2_triangulation_2<Gt, Tds>& p2t2, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P2T2, class GSOptions>
void add_to_graphics_scene(const P2T2& p2t2,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawPeriodic2Triangulation2

A shortcut to `CGAL::add_to_graphics_scene(p2t2, gs, Graphics_scene_options_periodic_2_triangulation_2{})`.
*/
template<class P2T2>
void add_to_graphics_scene(const P2T2& p2t2,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */
