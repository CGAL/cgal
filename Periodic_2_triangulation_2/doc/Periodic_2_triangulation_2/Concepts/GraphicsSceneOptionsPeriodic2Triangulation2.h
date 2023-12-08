/*!
\ingroup PkgPeriodic2Triangulation2Concepts

The concept `GraphicsSceneOptionsPeriodic2Triangulation2` defines data and methods used to tune the way that the cells of a `Periodic_2_triangulation_2` are considered for drawing or to be added into a graphics scene.

\cgalRefines{GraphicsSceneOptions}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Graphics_scene_options_periodic_2_triangulation_2 `CGAL::Graphics_scene_options_periodic_2_triangulation_2<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>`\endlink}
\cgalHasModelsEnd

*/
class GraphicsSceneOptionsPeriodic2Triangulation2
{
public:
  /// returns `true` if the domain of the Periodic_2_triangulation_2 must be drawn, `false` otherwise.
  /// Returns `false` by default.
  bool draw_domain() const;

  /// sets the draw domain value to `b`.
  void draw_domain(bool b);

  /// toggles the draw domain value.
  void toggle_draw_domain();

  /// returns the type of the display (STORED, UNIQUE, STORED_COVER_DOMAIN or UNIQUE_COVER_DOMAIN).
  typename DS::Iterator_type display_type() const;

  /// set the display type to the next type (in the ordered circular list STORED, UNIQUE, STORED_COVER_DOMAIN, UNIQUE_COVER_DOMAIN).
  void increase_display_type();

  /// returns the color used to draw the domain.
  const CGAL::IO::Color& domain_color() const;

  /// sets the color used to draw the domain to `c`.
  void domain_color(const CGAL::IO::Color& c)
};
