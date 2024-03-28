/*!
\ingroup PkgVoronoiDiagram2Concepts

The concept `GraphicsSceneOptionsVoronoiDiagram2` defines data and methods used to tune the way that the cells of a `Voronoi_diagram_2` are considered for drawing or to be added into a graphics scene.

\cgalRefines{GraphicsSceneOptions}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Graphics_scene_options_voronoi_diagram_2 `CGAL::Graphics_scene_options_voronoi_diagram_2<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>`\endlink}
\cgalHasModelsEnd

*/
class GraphicsSceneOptionsVoronoiDiagram2
{
public:
  /// returns the color of the dual vertices.
  const CGAL::IO::Color& dual_vertex_color() const;
  /// sets the color of dual vertices to `c`.
  void dual_vertex_color(const CGAL::IO::Color& c);

  /// returns the color of rays.
  const CGAL::IO::Color& ray_color() const;
  /// sets the color of rays to `c`.
  void ray_color(const CGAL::IO::Color& c);

  /// returns the color of the bisectors.
  const CGAL::IO::Color& bisector_color() const;
  /// sets the color of bisectors to `c`.
  void bisector_color(const CGAL::IO::Color& c);

  /// returns `true` if the voronoi vertices must be drawn, `false` otherwise.
  /// Returns `false` by default.
  bool draw_voronoi_vertices() const;
  /// sets the draw of voronoi vertices to `b`.
  void draw_voronoi_vertices(bool b);
  /// toggles the draw voronoi vertices value.
  void toggle_draw_voronoi_vertices();

  /// returns `true` if the dual vertices must be drawn, `false` otherwise.
  /// Returns `false` by default.
  bool draw_dual_vertices() const;
  /// sets the draw of dual vertices to `b`.
  void draw_dual_vertices();
  /// toggles the draw dual vertices value.
  void toggle_draw_dual_vertices();
};
