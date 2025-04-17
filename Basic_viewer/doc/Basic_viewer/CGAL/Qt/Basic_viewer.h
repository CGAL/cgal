namespace CGAL {
  namespace Qt {

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `Basic_viewer` is a Qt widget based on `QGLViewer` that allows to visualize 3D elements: points, segments, triangles, rays and lines. This class stores a reference to a `Graphics_scene`. Elements are added through the scene. This class requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER_QT` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer_Qt` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER_QT`.

CGAL::QGLViewer is our internal fork of <a href="https://github.com/GillesDebunne/libQGLViewer">QGLViewer class</a> which is <em>A versatile 3D OpenGL viewer based on QOpenGLWidget</em>.
*/
class Basic_viewer : public CGAL::QGLViewer
{
public:
  /// \name Constructors
  /// @{

  /// Constructor given a pointer on a `QWidget` (can be a `nullptr`) and a `Graphics_scene`.
  /// `title` will be the title of the window.
  Basic_viewer(QWidget* parent,
               const Graphics_scene& scene,
               const char* title="");

  /// @}

  /// \name Setters
  /// @{

  /// \brief Set size of vertices. 
  /// \param s The size of vertices.
  void size_vertices(float s);

  /// \brief Set size of edges. 
  /// \param s The size of edges.
  void size_edges(float s); 

  /// \brief Set size of rays. 
  /// \param s The size of rays.
  void size_rays(float s); 

  /// \brief Set size of lines. 
  /// \param s The size of lines.
  void size_lines(float s); 

  /// \brief Enables or disables the drawing of vertices.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_vertices(bool b);

  /// \brief Enables or disables the drawing of edges.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_edges(bool b);

  /// \brief Enables or disables the drawing of rays.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_rays(bool b);

  /// \brief Enables or disables the drawing of lines.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_lines(bool b);

  /// \brief enables or disables the drawing of faces.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_faces(bool b);

  /// enables or disables the drawing of texts.
  /// \brief enables or disables the drawing of texts.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_text(bool b);

  /// \brief Enables or disables the drawing of mesh triangles.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_mesh_triangles(bool b);

  /// \brief Enables or disables the use of only one color or the use of multiple colors.
  /// \param b Set to `true` to use only one color, `false` to use multiple colors.
  void use_default_color(bool b);

  /// \brief Enables or disables the use of a single color for all normals.
  /// \param b Set to `true` to enable, `false` to disable.
  void use_default_color_normals(bool b);

  /// \brief enables or disables the use of flat shading or the use of smooth shading.
  /// \param b Set to `true` for flat shading, `false` for smooth shading.
  void flat_shading(bool b);

  /// \brief Enables or disables the reversal of normals.
  /// \param b Set to `true` to reverse normals, `false` to keep normals as is.
  void reverse_normal(bool b);

  /// \brief Sets the default color of the normals.
  /// \param c The default color of the normals.
  void default_color_normals(const CGAL::IO::Color& c);

  /// \brief Sets the height factor value of the normals.
  /// \param h The height factor value of the normals.
  void normal_height_factor(float h);

  /// \brief Toggles the drawing of vertices.
  void toggle_draw_vertices();

  /// \brief Toggles the drawing of edges.
  void toggle_draw_edges();

  /// \brief Toggles the drawing of rays.
  void toggle_draw_rays();

  /// \brief Toggles the drawing of lines.
  void toggle_draw_lines();

  /// \brief Toggles the drawing of faces.
  void toggle_draw_faces();

  /// \brief Toggles the use of mono color mode.
  void toggle_use_default_color();

  /// \brief Toggles the use of the default color mode for normals.
  void toggle_use_default_color_normal();

  /// \brief Toggles the use of flat shading.
  void toggle_flat_shading();

  /// \brief Toggles the drawing of text.
  void toggle_draw_text();

  /// \brief Reverses all normals of vertices and faces.
  void reverse_all_normals();

  /// @}

  /// \name Getters
  /// @{

  /// \brief Checks if vertices are drawn.
  /// \return `true` if vertices are drawn, `false` otherwise.
  bool draw_vertices() const;

  /// \brief Checks if edges are drawn.
  /// \return `true` if edges are drawn, `false` otherwise.
  bool draw_edges() const;

  /// \brief Checks if rays are drawn.
  /// \return `true` if rays are drawn, `false` otherwise.
  bool draw_rays() const;

  /// \brief Checks if lines are drawn.
  /// \return `true` if lines are drawn, `false` otherwise.
  bool draw_lines() const;

  /// \brief Checks if faces are drawn.
  /// \return `true` if faces are drawn, `false` otherwise.
  bool draw_faces() const;

  /// \brief Checks if text is drawn.
  /// \return `true` if text is drawn, `false` otherwise.
  bool draw_text() const;

  /// \brief Checks if the default color mode is used.
  /// \return `true` if mono color mode is used, `false` otherwise.
  bool use_default_color() const;

  /// \brief Checks if the default color mode for normals is used.
  /// \return `true` if default color mode for normals is used, `false` otherwise.
  bool use_default_color_normal() const;

  /// \brief Checks if normals are reversed.
  /// \return `true` if normals are reversed, `false` otherwise.
  bool reverse_normal() const;

  /// \brief Checks if the clipping plane is enabled.
  /// \return `true` if the clipping plane is enabled, `false` otherwise.
  bool clipping_plane_enabled() const;

  /// \brief Checks if m_no_2D_mode is false and the graphics scene is two-dimensional.
  /// \return `true` if m_no_2D_mode is false and the scene is 2D, `false` otherwise.
  bool is_two_dimensional() const;

  /// \brief Gets the clipping plane when enabled.
  /// \return The clipping plane as a `CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3` object.
  CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 clipping_plane() const;

  /// \brief Gets the graphics scene of the viewer.
  /// \return A reference to the `Graphics_scene` object.
  const Graphics_scene& graphics_scene() const;

  /// @}

  /// \name Draw
  /// @{

  /// draws the viewer without recomputing all internal buffers.
  virtual void draw();

  /// redraws the viewer, i.e., recompute all internal buffers and update the window.
  virtual void redraw();

  /// @}

  /// Function called when a key is pressed. Users can define their own function in order
  /// to add specific behavior.
  std::function<bool(QKeyEvent *, CGAL::Qt::Basic_viewer *)> on_key_pressed;
};


//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `QApplication_and_basic_viewer` regroups a `Basic_viewer` and Qt `QApplication`. The `QApplication` is created in the constructor, but started by the `run()` method. This allows for example users to modify the `on_key_pressed` method of the `Basic_viewer` to define their own behavior.

This class requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER_QT` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer_Qt` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER_QT`.

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

 opens a new window and draws the given `Graphics_scene` (which must have been filled before). `title` will be the title of the window.  A call to this method is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER_QT` is defined. Linking with the cmake target `CGAL::CGAL_Basic_viewer_Qt` will link with `CGAL_Qt6` and add the definition CGAL_USE_BASIC_VIEWER_QT.
*/
void draw_graphics_scene(const Graphics_scene& graphic_scene,
                         const char *title="CGAL Basic Viewer")
{}

} // End namespace CGAL

