namespace CGAL {

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `Basic_viewer_qt` is a Qt widget based on `QGLViewer` that allows to visualize 3D elements: points, segments, triangles, rays and lines. This class stores a reference to a `Graphics_scene`. Elements can either be added directly in the viewer or through the scene. This class requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam BufferType the number type used for point coordinates: `float` by default.

*/
template <typename BufferType=float>
class Basic_viewer_qt : public CGAL::QGLViewer
{
public:
  /// Constructor given a pointer on a `QWidget` (can be a `nullptr`) and a `Graphics_scene`.
  /// `title` will be the title of the window.
  Basic_viewer_qt(QWidget* parent,
                  Graphics_scene<BufferType>& buf,
                  const char* title="");

  /// enables or disables the drawing of vertices.
  void set_draw_vertices(bool b);

  /// enables or disables the drawing of edges.
  void set_draw_edges(bool b);

  /// enables or disables the drawing of rays.
  void set_draw_rays(bool b);

  /// enables or disables the drawing of lines.
  void set_draw_lines(bool b);

  /// enables or disables the drawing of faces.
  void set_draw_faces(bool b);

  /// enables or disables the use of only one color (if `b` is `true`) or the use of multiple colors (if `b` is `false`).
  void set_use_mono_color(bool b);

  /// enables or disables the drawing of texts.
  void set_draw_text(bool b);

  /// sets the color used for vertices in mono_color mode.
  void set_vertices_mono_color(const CGAL::IO::Color& c);

  /// sets the color used for edges in mono_color mode.
  void set_edges_mono_color(const CGAL::IO::Color& c);

  /// sets the color used for rays in mono_color mode.
  void set_rays_mono_color(const CGAL::IO::Color& c);

  /// sets the color used for lines in mono_color mode.
  void set_lines_mono_color(const CGAL::IO::Color& c);

  /// sets the color used for faces in mono_color mode.
  void set_faces_mono_color(const CGAL::IO::Color& c);

  /// negates the drawing of vertices (becomes `true` if it was `false` and reciprocally).
  void negate_draw_vertices();

  /// negates the drawing of edges (becomes `true` if it was `false` and reciprocally).
  void negate_draw_edges();

  /// negates the drawing of rays (becomes `true` if it was `false` and reciprocally).
  void negate_draw_rays();

  /// negates the drawing of lines (becomes `true` if it was `false` and reciprocally).
  void negate_draw_lines();

  /// negates the drawing of faces (becomes `true` if it was `false` and reciprocally).
  void negate_draw_faces();

  /// negates the use of mono color mode (becomes `true` if it was `false` and reciprocally).
  void negate_use_mono_color();

  /// negates the drawing of text (becomes `true` if it was `false` and reciprocally).
  void negate_draw_text();

  /// returns `true` if vertices are drawn.
  bool get_draw_vertices() const;

  /// returns `true` if edges are drawn.
  bool get_draw_edges() const;

  /// returns `true` if rays are drawn.
  bool get_draw_rays() const;

  /// returns `true` if lines are drawn.
  bool get_draw_lines() const;

  /// returns `true` if faces are drawn.
  bool get_draw_faces() const;

  /// returns `true` if mono color mode is used.
  bool get_use_mono_color() const;

  /// returns `true` if normal are reversed.
  bool get_inverse_normal() const;

  /// returns `true` if text are drawn.
  bool get_draw_text() const;

  /// returns the mono color used for vertices.
  const CGAL::IO::Color& get_vertices_mono_color() const;

  /// returns the mono color used for edges.
  const CGAL::IO::Color& get_edges_mono_color() const;

  /// returns the mono color used for rays.
  const CGAL::IO::Color& get_rays_mono_color() const;

  /// returns the mono color used for lines.
  const CGAL::IO::Color& get_lines_mono_color() const;

  /// returns the mono color used for faces.
  const CGAL::IO::Color& get_faces_mono_color() const;

  /// clears the basic viewer, i.e., remove all its elements.
  void clear();

  /// returns `true` if the viewer is empty.
  bool is_empty() const;

  /// returns the bounding box of all the elements in the viewer.
  const CGAL::Bbox_3& bounding_box() const;

  /// returns `true` if the clipping plane is enabled.
  bool is_clipping_plane_enabled() const;

  /// returns the clipping plane when it is enabled.
  CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 clipping_plane() const;

  /// adds the given point in the viewer.
  template <typename KPoint>
  void add_point(const KPoint &p);

  /// adds the given colored point in the viewer.
  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor);

  /// adds the given segment in the viewer.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2);

  /// adds the given colored segment in the viewer.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor);

  /// adds the given ray in the viewer: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v);

  /// adds the given colored ray in the viewer: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor);

  /// adds the given line in the viewer, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v);

  /// adds the given colored line in the viewer, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor);

  /// starts a new face.
  void face_begin();

  /// starts a new colored face.
  void face_begin(const CGAL::IO::Color &acolor);

  /// returns `true` iff a face is started.
  bool is_a_face_started() const;

  /// adds the given point in the current face.
  /// @pre `is_a_face_started()`
  template <typename KPoint> bool add_point_in_face(const KPoint &kp);

  /// adds the given point in the current face, having the vertex normal.
  /// @pre `is_a_face_started()`
  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal);

  /// ends the current face.
  /// @pre `is_a_face_started()`
  void face_end();

  /// adds the given text at the given position in the viewer.
  template <typename KPoint>
  void add_text(const KPoint &kp, const char *txt);

  /// adds the given text at the given position in the viewer.
  template <typename KPoint>
  void add_text(const KPoint &kp, const std::string &txt);

  /// returns the graphics scene of the viewer.
  Graphics_scene<BufferType>& get_graphics_scene();

  /// returns the graphics scene of the viewer, const version.
  const Graphics_scene<BufferType>& get_graphics_scene() const;

  /// negates all normals of vertices and faces.
  void negate_all_normals();

  /// returns `true` if the data structure in in 2D, i.e., lies on a plane.
  bool is_two_dimensional() const;

  /// draws the viewer without recomputing all internal buffers.
  virtual void draw();

  /// redraws the viewer, i.e., recompute all internal buffers and update the window.
  virtual void redraw();

  /// Function called when a key is pressed. Users can define their own function in order
  /// to add specific behavior.
  std::function<bool(QKeyEvent *, CGAL::Basic_viewer_qt<BufferType> *)> on_key_pressed;
};


//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `QApplication_and_basic_viewer` regroups a `Basic_viewer_qt` and Qt `QApplication`. The `QApplication` is created in the constructor, but started by the `run()` method. This allows for example users to modify the `on_key_pressed` method of the `Basic_viewer_qt` to define their own behavior. This class requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam BufferType the number type used for point coordinates: `float` by default.

*/
template <typename BufferType=float>
class QApplication_and_basic_viewer
{
public:
  /// Constructor given a `Graphics_scene` and possibly a title.
  QApplication_and_basic_viewer(CGAL::Graphics_scene<BufferType>& buffer,
                                const char* title="CGAL Basic Viewer");

  /// runs the `QApplication`, i.e., open the Qt window. A call to this method is blocking, that is the program continues as soon as the user closes the window.
  void run();

  /// returns a reference to the `Basic_viewer_qt` associated with this.
  Basic_viewer_qt<BufferType>& basic_viewer();
};

//------------------------------------------------------------------------------
/*!
  opens a new window and draws the given `Graphics_scene` (which must have been filled before). `title` will be the title of the window.  A call to this method is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition CGAL_USE_BASIC_VIEWER.
*/
template <typename BufferType=float>
void draw_graphics_scene(Graphics_scene<BufferType>& graphic_buffer,
                         const char *title="CGAL Basic Viewer")
{}

} // End namespace CGAL
