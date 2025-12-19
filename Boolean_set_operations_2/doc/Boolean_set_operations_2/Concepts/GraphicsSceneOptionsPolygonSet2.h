/*!
\ingroup PkgBooleanSetOperations2Concepts

The concept `GraphicsSceneOptionsPolygonSet2` defines data and methods used to tune the way that the cells of a `Polygon_set_2` are considered for drawing or to be added into a graphics scene.

\cgalRefines{GraphicsSceneOptions}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Graphics_scene_options_polygon_set_2 `CGAL::Graphics_scene_options_polygon_set_2<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>`\endlink}
\cgalHasModelsEnd

*/
class GraphicsSceneOptionsPolygonSet2
{
public:
  /// returns the color of the unbounded face.
  const CGAL::IO::Color& unbounded_face_color() const;
  /// sets the color of the unbounded face to `c`.
  void unbounded_face_color(const CGAL::IO::Color& c);

  /// returns `true` if the unbounded face must be drawn, `false` otherwise.
  /// Returns `false` by default.
  bool draw_unbounded() const;
  /// sets the draw of unbounded face to `b`.
  void draw_unbounded(bool b);
  /// toggles the draw unbounded face value.
  void toggle_draw_unbounded();

  /// returns the height of the box used to draw the unbounded face.
  int height() const;
  /// returns the width of the box used to draw the unbounded face.
  int width() const;

  /// sets the height of the box used to draw the unbounded face to `i`.
  void height(int i);
  /// sets the width of the box used to draw the unbounded face to `i`.
  void width(int i);
};
