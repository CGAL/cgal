namespace CGAL {

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `Basic_viewer_qt` is a Qt widget based on `QGLViewer` that allows to visualize 3D elements: points, segments, triangles, rays and lines. This class stores a reference to a `Graphic_storage`. Elements can either be added directly in the viewer or through the storage. This class requires CGAL_Qt5, and is only available if the macro CGAL_USE_BASIC_VIEWER is defined. Linking with the cmake target CGAL::CGAL_Basic_viewer will link with CGAL_Qt5 and add the definition CGAL_USE_BASIC_VIEWER.

\tparam BufferType the type used for point coordinates: `float` by default.

*/
template <typename BufferType=float>
class Basic_viewer_qt : public CGAL::QGLViewer
{
public:
  /// constructor given a pointer on a `QWidget` (can be a `nullptr`) and a `Graphic_storage`.
  /// `title` will be the title of the window.
  Basic_viewer_qt(QWidget* parent,
                  Graphic_storage<BufferType>& buf,
                  const char* title="");

  /// enable or disable the drawing of vertices.
  void set_draw_vertices(bool b);

  /// enable or disable the drawing of edges.
  void set_draw_edges(bool b);

  /// enable or disable the drawing of rays.
  void set_draw_rays(bool b);

  /// enable or disable the drawing of lines.
  void set_draw_lines(bool b);

  /// enable or disable the drawing of faces.
  void set_draw_faces(bool b);

  /// enable or disable the use of only one color (if `b` is `true`) or the use of multiple colors (if `b` is `false`).
  void set_use_mono_color(bool b);

  /// enable or disable the drawing of texts.
  void set_draw_text(bool b);

  /// set the color used for vertices in mono_color mode.
  void set_vertices_mono_color(const CGAL::IO::Color& c);

  /// set the color used for edges in mono_color mode.
  void set_edges_mono_color(const CGAL::IO::Color& c);

  /// set the color used for rays in mono_color mode.
  void set_rays_mono_color(const CGAL::IO::Color& c);

  /// set the color used for lines in mono_color mode.
  void set_lines_mono_color(const CGAL::IO::Color& c);

  /// set the color used for faces in mono_color mode.
  void set_faces_mono_color(const CGAL::IO::Color& c);

  /// negate the drawing of vertices (becomes `true` if it was `false` and reciprocally).
  void negate_draw_vertices();

  /// negate the drawing of edges (becomes `true` if it was `false` and reciprocally).
  void negate_draw_edges();

  /// negate the drawing of rays (becomes `true` if it was `false` and reciprocally).
  void negate_draw_rays();

  /// negate the drawing of lines (becomes `true` if it was `false` and reciprocally).
  void negate_draw_lines();

  /// negate the drawing of faces (becomes `true` if it was `false` and reciprocally).
  void negate_draw_faces();

  /// negate the use of mono color mode (becomes `true` if it was `false` and reciprocally).
  void negate_use_mono_color();

  /// negate the drawing of text (becomes `true` if it was `false` and reciprocally).
  void negate_draw_text();

  /// @return `true` if vertices are drawn.
  bool get_draw_vertices() const;

  /// @return `true` if edges are drawn.
  bool get_draw_edges() const;

  /// @return `true` if rays are drawn.
  bool get_draw_rays() const;

  /// @return `true` if lines are drawn.
  bool get_draw_lines() const;

  /// @return `true` if faces are drawn.
  bool get_draw_faces() const;

  /// @return `true` if mono color mode is used.
  bool get_use_mono_color() const;

  /// @return `true` if normal are reversed.
  bool get_inverse_normal() const;

  /// @return `true` if text are drawn.
  bool get_draw_text() const;

  /// @return the mono color used for vertices.
  const CGAL::IO::Color& get_vertices_mono_color() const;

  /// @return the mono color used for edges.
  const CGAL::IO::Color& get_edges_mono_color() const;

  /// @return the mono color used for rays.
  const CGAL::IO::Color& get_rays_mono_color() const;

  /// @return the mono color used for lines.
  const CGAL::IO::Color& get_lines_mono_color() const;

  /// @return the mono color used for faces.
  const CGAL::IO::Color& get_faces_mono_color() const;

  /// clear the basic viewer, i.e. remove all its elements.
  void clear();

  /// @return `true` if the viewer is empty.
  bool is_empty() const;

  /// @return the bounding box of all the elements in the viewer.
  const CGAL::Bbox_3& bounding_box() const;

  /// @return `true` if the clipping plane is enabled.
  bool is_clipping_plane_enabled() const;

  /// @return the clipping plane when it is enabled.
  Local_kernel::Plane_3 clipping_plane() const;

  /// add the given point in the viewer.
  template <typename KPoint>
  void add_point(const KPoint &p);

  /// add the given colored point in the viewer.
  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor);

  /// add the given segment in the viewer.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2);

  /// add the given colored segment in the viewer.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor);

  /// add the given ray in the viewer: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v);

  /// add the given colored ray in the viewer: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor);

  /// add the given line in the viewer, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v);

  /// add the given colored line in the viewer, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor);

  /// start a new face.
  void face_begin();

  /// start a new colored face.
  void face_begin(const CGAL::IO::Color &acolor);

  /// @return `true` iff a face is started.
  bool is_a_face_started() const;

  /// add the given point in the current face.
 /// @pre `is_a_face_started()`
  template <typename KPoint> bool add_point_in_face(const KPoint &kp);

  /// add the given point in the current face, having the vertex normal.
 /// @pre `is_a_face_started()`
  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal);

  /// end the current face.
 /// @pre `is_a_face_started()`
  void face_end();

  /// add the given text at the given position in the viewer.
  template <typename KPoint>
  void add_text(const KPoint &kp, const char *txt);

  /// add the given text at the given position in the viewer.
  template <typename KPoint>
  void add_text(const KPoint &kp, const std::string &txt);

  /// @return the graphic storage of the viewer.
  Graphic_storage<BufferType>& get_graphic_storage();

  /// @return the graphic storage of the viewer, const version.
  const Graphic_storage<BufferType>& get_graphic_storage() const;

  /// negate all normal of vertices and faces.
  void negate_all_normals();

  // \return true if the data structure in in 2D, i.e. lies on a plane.
  bool is_two_dimensional() const;

  /// draw the viewer without recomputing all internal buffers.
  virtual void draw();

  /// redraw the viewer, i.e. recompute all internal buffers and update the window.
  virtual void redraw();

  /// function called when a key is pressed. Users can define their own function in order
  /// to add specific behavior.
  std::function<bool(QKeyEvent *, CGAL::Basic_viewer_qt<BufferType> *)> on_key_pressed;
};


//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `QApplication_and_basic_viewer` regroups a `Basic_viewer_qt` and Qt `QApplication`. The `QApplication` is created in the constructor, but ran by the run method. This allows for example users to modify the `on_key_pressed` method of the `Basic_viewer_qt` to define their own behavior. This class requires CGAL_Qt5, and is only available if the macro CGAL_USE_BASIC_VIEWER is defined. Linking with the cmake target CGAL::CGAL_Basic_viewer will link with CGAL_Qt5 and add the definition CGAL_USE_BASIC_VIEWER.

\tparam BufferType the type used for point coordinates: `float` by default.

*/
template <typename BufferType=float>
class QApplication_and_basic_viewer
{
public:
  /// Constructor given a Graphic_storage and possibly a title.
  QApplication_and_basic_viewer(CGAL::Graphic_storage<BufferType>& buffer,
                                const char* title="CGAL Basic Viewer");

  /// run the QApplication, i.e. open the Qt window. A call to this method is blocking, that is the program continues as soon as the user closes the window.
  void run();

  /// @return a reference to the `Basic_viewer_qt` associated with this.
  Basic_viewer_qt<BufferType>& basic_viewer();
};

//------------------------------------------------------------------------------
/*!
  opens a new window and draws the given `Graphic_storage` (which must have been filled before). `title` will be the title of the window.  A call to this method is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the macro CGAL_USE_BASIC_VIEWER is defined. Linking with the cmake target CGAL::CGAL_Basic_viewer will link with CGAL_Qt5 and add the definition CGAL_USE_BASIC_VIEWER.
*/
template <typename BufferType=float>
void draw_graphic_storage(Graphic_storage<BufferType>& graphic_buffer,
                          const char *title="CGAL Basic Viewer")
{}

} // End namespace CGAL

