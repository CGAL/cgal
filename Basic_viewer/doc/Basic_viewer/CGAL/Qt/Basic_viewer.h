namespace CGAL {
  namespace Qt {

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `Basic_viewer` is a Qt widget based on `QGLViewer` that allows to visualize 3D elements: points, segments, triangles, rays and lines. This class stores a reference to a `Graphics_scene`. Elements are added through the scene. This class requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER_QT` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` or `CGAL::CGAL_Basic_viewer_Qt` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER` & `CGAL_USE_BASIC_VIEWER_QT`.

CGAL::QGLViewer is our internal fork of <a href="https://github.com/GillesDebunne/libQGLViewer">QGLViewer class</a> which is <em>A versatile 3D OpenGL viewer based on QOpenGLWidget</em>.
*/
class Basic_viewer : public CGAL::QGLViewer
{
public:
  /// Constructor given a pointer on a `QWidget` (can be a `nullptr`) and a `Graphics_scene`.
  /// `title` will be the title of the window.
  Basic_viewer(QWidget* parent,
               const Graphics_scene& scene,
               const char* title="");

  /// Set size of points 
  void size_points(float s);

  /// Set size of edges 
  void size_edges(float s); 

  /// Set size of rays 
  void size_rays(float s); 

  /// Set size of lines 
  void size_lines(float s); 

  /// enables or disables the drawing of vertices.
  void draw_vertices(bool b);

  /// enables or disables the drawing of edges.
  void draw_edges(bool b);

  /// enables or disables the drawing of rays.
  void draw_rays(bool b);

  /// enables or disables the drawing of lines.
  void draw_lines(bool b);

  /// enables or disables the drawing of faces.
  void draw_faces(bool b);

  /// enables or disables the use of only one color (if `b` is `true`) or the use of multiple colors (if `b` is `false`).
  void use_default_color(bool b);

  /// enables or disables the use of flat shading (if `b` is `true`) or the use of smooth shading (if `b` is `false`).
  void flat_shading(bool b);

  /// enables or disables the drawing of texts.
  void draw_text(bool b);

  /// toggles the drawing of vertices.
  void toggle_draw_vertices();

  /// toggles the drawing of edges.
  void toggle_draw_edges();

  /// toggles the drawing of rays.
  void toggle_draw_rays();

  /// toggles the drawing of lines.
  void toggle_draw_lines();

  /// toggles the drawing of faces.
  void toggle_draw_faces();

  /// toggles the use of mono color mode.
  void toggle_use_default_color();

  /// toggles the use of flat shading.
  void toggle_flat_shading();

  /// toggles the drawing of text.
  void toggle_draw_text();

  /// returns `true` if vertices are drawn.
  bool draw_vertices() const;

  /// returns `true` if edges are drawn.
  bool draw_edges() const;

  /// returns `true` if rays are drawn.
  bool draw_rays() const;

  /// returns `true` if lines are drawn.
  bool draw_lines() const;

  /// returns `true` if faces are drawn.
  bool draw_faces() const;

  /// returns `true` if mono color mode is used.
  bool use_default_color() const;

  /// returns `true` if normals are reversed.
  bool reverse_normal() const;

  /// returns `true` if text are drawn.
  bool draw_text() const;

  /// returns `true` if the clipping plane is enabled.
  bool clipping_plane_enabled() const;

  /// returns the clipping plane when it is enabled.
  CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 clipping_plane() const;

  /// returns the graphics scene of the viewer.
  const Graphics_scene& graphics_scene() const;

  /// reverses all normals of vertices and faces.
  void reverse_all_normals();

  /// draws the viewer without recomputing all internal buffers.
  virtual void draw();

  /// redraws the viewer, i.e., recompute all internal buffers and update the window.
  virtual void redraw();

  /// Function called when a key is pressed. Users can define their own function in order
  /// to add specific behavior.
  std::function<bool(QKeyEvent *, CGAL::Qt::Basic_viewer *)> on_key_pressed;
};


//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `QApplication_and_basic_viewer` regroups a `Basic_viewer` and Qt `QApplication`. The `QApplication` is created in the constructor, but started by the `run()` method. This allows for example users to modify the `on_key_pressed` method of the `Basic_viewer` to define their own behavior.

This class requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` or `CGAL_USE_BASIC_VIEWER_QT` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` or `CGAL::CGAL_Basic_viewer_Qt` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER` & `CGAL_USE_BASIC_VIEWER_QT`.

*/
class QApplication_and_basic_viewer
{
public:
  /// Constructor given a `Graphics_scene` and possibly a title.
  QApplication_and_basic_viewer(const CGAL::Graphics_scene& gs,
                                const char* title="CGAL Basic Viewer");

  /// runs the `QApplication`, i.e., open the Qt window. A call to this method is blocking, that is the program continues as soon as the user closes the window.
  void run();

  /// returns a reference to the `Basic_viewer` associated with this.
  Basic_viewer& basic_viewer();
};

} // End namespace Qt

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

 opens a new window and draws the given `Graphics_scene` (which must have been filled before). `title` will be the title of the window.  A call to this method is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition CGAL_USE_BASIC_VIEWER.
*/
void draw_graphics_scene(const Graphics_scene& graphic_scene,
                         const char *title="CGAL Basic Viewer")
{}

} // End namespace CGAL

