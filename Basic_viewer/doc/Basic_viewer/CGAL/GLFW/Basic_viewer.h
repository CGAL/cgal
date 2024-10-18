namespace CGAL {
  namespace GLFW {

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

The class `Basic_viewer` is an OpenGL (GLFW) viewer that allows visualizing 3D elements such as points, segments, triangles, rays, and lines. This class stores a reference to a `Graphics_scene`, which contains the elements to be visualized. These elements are added through the scene. The class requires `CGAL_GLFW` and is only available if the macro `CGAL_USE_BASIC_VIEWER_GLFW` is defined. Linking with the CMake target `CGAL::CGAL_Basic_viewer_GLFW` will automatically link with `CGAL_GLFW` and add the definition `CGAL_USE_BASIC_VIEWER_GLFW`.

This viewer is designed to be lightweight and easy to use, providing essential 3D visualization capabilities for developers working with CGAL. The input is managed by a class based on GLFW, allowing for cross-platform interaction.
*/
class Basic_viewer : public Input
{
public:
  /// \name Constructors
  /// @{

  /*!
    \brief Constructs a `Basic_viewer` given a `Graphics_scene` and an optional window title.

    \param scene The graphics scene containing the 3D elements to visualize.
    \param title The title of the viewer window. Defaults to an empty string.
  */
  Basic_viewer(const Graphics_scene& scene,
               const char* title="");

  /// @}

  /// \name Initialization
  /// @{

  /*!
    \brief Initializes the OpenGL context, window, callbacks, and buffers.

    This method is called in the constructor and sets up everything required for the viewer to function. If you intend to use `make_screenshot()` instead of `show()`, you can pass `screenshotOnly=true` to avoid unnecessary initializations such as the event handling system.

    \param screenshotOnly If set to `true`, skips certain initializations. Defaults to `false`.
  */
  void initialize(bool screenshotOnly=false);

  /// @}

  /// \name Main Application Loop
  /// @{

  /*!
    \brief Starts the application loop to draw the scene and handle events.

    This method blocks the execution, keeping the window open until the user closes it. During this loop, the viewer continuously renders the scene and responds to user input.
  */
  void show();

  /*!
    \brief Draws a single frame and saves it as an image file.

    This function allows you to capture the current view of the scene and save it to a specified file. It can be used to generate screenshots without displaying the interactive viewer.

    \param filePath The path where the screenshot will be saved.
  */
  void make_screenshot(const std::string& filePath);

  /// @}

  /// \name Setters
  /// @{

  /// \brief Sets the window size.
  /// \param size The size of the window as a 2D vector.
  void window_size(const Eigen::Vector2f &s);

  /// \brief Sets the size of vertices.
  /// \param s The size of vertices.
  void size_vertices(float s);

  /// \brief Sets the size of edges.
  /// \param s The size of edges.
  void size_edges(float s); 

  /// \brief Sets the size of rays.
  /// \param s The size of rays.
  void size_rays(float s); 

  /// \brief Sets the size of lines.
  /// \param s The size of lines.
  void size_lines(float s);

  /// \brief Sets the position of the light.
  /// \param p The position of the light as a 4D vector.
  void light_position(const Eigen::Vector4f& p);

  /// \brief Sets the ambient color of the light.
  /// \param c The ambient color as a 4D vector.
  void light_ambient(const Eigen::Vector4f& c);

  /// \brief Sets the diffuse color of the light.
  /// \param c The diffuse color as a 4D vector.
  void light_diffuse(const Eigen::Vector4f& c);

  /// \brief Sets the specular color of the light.
  /// \param c The specular color as a 4D vector.
  void light_specular(const Eigen::Vector4f& c);

  /// \brief Sets the shininess/brightness value of the light.
  /// \param s The shininess value.
  void light_shininess(float s);
  
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

  /// \brief Enables or disables the drawing of faces.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_faces(bool b);

  /// \brief Enables or disables the clipping plane.
  /// \param b Set to `true` to enable the clipping plane, `false` to disable.
  void draw_clipping_plane(bool b);

  /// \brief Enables or disables the drawing of the world axis.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_world_axis(bool b);

  /// \brief Enables or disables the drawing of the XY grid.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_xy_grid(bool b);

  /// \brief Enables or disables the drawing of mesh triangles.
  /// \param b Set to `true` to enable, `false` to disable.
  void draw_mesh_triangles(bool b);

  /// \brief Enables or disables the use of a single color for all elements.
  /// \param b Set to `true` to use a single color, `false` for multiple colors.
  void use_default_color(bool b);

  /// \brief Enables or disables the use of a single color for all normals.
  /// \param b Set to `true` to enable, `false` to disable.
  void use_default_color_normals(bool b);

  /// \brief Enables or disables the printing of the application state to the terminal.
  /// \param b Set to `true` to enable printing, `false` to disable.
  void print_application_state(bool b);

  /// \brief Enables or disables the reversal of normals.
  /// \param b Set to `true` to reverse normals, `false` to keep normals as is.
  void reverse_normal(bool b);

  /// \brief Enables or disables flat shading (as opposed to smooth shading).
  /// \param b Set to `true` for flat shading, `false` for smooth shading.
  void flat_shading(bool b);

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

  /// \brief Toggles the drawing of the clipping plane.
  void toggle_draw_clipping_plane();

  /// \brief Toggles the drawing of the world axis.
  void toggle_draw_world_axis();

  /// \brief Toggles the drawing of the XY grid.
  void toggle_draw_xy_grid(); 

  /// \brief Toggles the drawing of mesh triangles.
  void toggle_draw_mesh_triangles(); 

  /// \brief Toggles the use of the default color mode.
  void toggle_use_default_color();

  /// \brief Toggles the use of the default color mode for normals.
  void toggle_use_default_color_normal();

  /// \brief Toggles the printing of the application state.
  void toggle_print_application_state();
  
  /// \brief Sets the clipping plane display mode.
  /// 
  /// \param m The display mode:
  /// - `CLIPPING_PLANE_OFF`: Disables the clipping plane.
  /// - `CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF`: Renders half of the clipped object as solid and half as transparent.
  /// - `CLIPPING_PLANE_SOLID_HALF_WIRE_HALF`: Renders half of the clipped object as solid and half with wireframe.
  /// - `CLIPPING_PLANE_SOLID_HALF_ONLY`: Only renders the positive side of the clipped object.
  void display_mode(DisplayMode m);

  /// \brief Enables the 2D/orthographic camera mode.
  void two_dimensional(); 

  /// \brief Sets the input manager to AZERTY keyboard layout.
  void azerty_layout();

  /// \brief Sets the duration of animations.
  /// \param d The duration of animations as a `std::chrono::milliseconds` object.
  void animation_duration(std::chrono::milliseconds d);

  /// \brief Sets the radius of the scene.
  /// \param r The radius of the scene.
  void scene_radius(float r);

  /// \brief Sets the center of the scene.
  /// \param x The X-coordinate of the center.
  /// \param y The Y-coordinate of the center.
  /// \param z The Z-coordinate of the center.
  void scene_center(float x, float y, float z);

  /// \brief Sets the center of the scene using a 3D vector.
  /// \param c The center of the scene as an Eigen vector.
  void scene_center(const Eigen::Vector3f& c);

  /// \brief Sets the position of the camera.
  /// \param x The X-coordinate of the camera position.
  /// \param y The Y-coordinate of the camera position.
  /// \param z The Z-coordinate of the camera position.
  void camera_position(float x, float y, float z);

  /// \brief Sets the position of the camera using a 3D vector.
  /// \param p The camera position as an Eigen vector.
  void camera_position(const Eigen::Vector3f& p);

  /// \brief Sets the orientation of the camera.
  /// \param f The forward direction of the camera as an Eigen vector.
  /// \param u The up vector angle (in degrees) of the camera.
  void camera_orientation(const Eigen::Vector3f& f, float u);

  /// \brief Sets the zoom level of the camera.
  /// \param z The zoom level.
  void zoom(float z);

  /// \brief Aligns the camera forward direction vector to the clipping plane normal.
  void align_camera_to_clipping_plane();

  /// \brief Sets the orientation of the clipping plane.
  /// \param n The normal vector of the clipping plane.
  void clipping_plane_orientation(const Eigen::Vector3f& n);

  /// \brief Translates the clipping plane along its normal.
  /// \param t The translation amount.
  void clipping_plane_translate_along_normal(float t);

  /// \brief Translates the clipping plane along the camera forward direction.
  /// \param t The translation amount.
  void clipping_plane_translate_along_camera_forward(float t);

  /// \brief Reverses all normals of vertices and faces.
  void reverse_all_normals();

  /// @}

  /// \name Getters
  /// @{

  /// \brief Gets the camera position.
  /// \return The camera position as an Eigen vector.
  Eigen::Vector3f position() const;

  /// \brief Gets the camera forward direction.
  /// \return The forward direction as an Eigen vector.
  Eigen::Vector3f forward() const;

  /// \brief Gets the camera right direction.
  /// \return The right direction as an Eigen vector.
  Eigen::Vector3f right() const;

  /// \brief Gets the camera up direction.
  /// \return The up direction as an Eigen vector.
  Eigen::Vector3f up() const;

  /// \brief Gets the size of points.
  /// \return The size of points.
  float size_vertices() const;

  /// \brief Gets the size of edges.
  /// \return The size of edges.
  float size_edges() const;
  
  /// \brief Gets the size of rays.
  /// \return The size of rays.
  float size_rays() const;

  /// \brief Gets the size of lines.
  /// \return The size of lines.
  float size_lines() const;

  /// \brief Gets the position of the light.
  /// \return The position of the light as a 4D vector.
  vec4f light_position() const;

  /// \brief Gets the ambient color of the light.
  /// \return The ambient color as a 4D vector.
  vec4f light_ambient() const;

  /// \brief Gets the diffuse color of the light.
  /// \return The diffuse color as a 4D vector.
  vec4f light_diffuse() const;

  /// \brief Gets the specular color of the light.
  /// \return The specular color as a 4D vector.
  vec4f light_specular() const;

  /// \brief Gets the shininess/brightness value of the light.
  /// \return The shininess value.
  float light_shininess() const;

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

  /// \brief Checks if world axis are drawn.
  /// \return `true` if world axis are drawn, `false` otherwise.
  bool draw_world_axis() const;

  /// \brief Checks if XY grid is drawn.
  /// \return `true` if XY grid is drawn, `false` otherwise.
  bool draw_xy_grid() const;

  /// \brief Checks if mesh triangles are drawn.
  /// \return `true` if mesh triangles are drawn, `false` otherwise.
  bool draw_mesh_triangles() const;

  /// \brief Checks if the default color mode is used.
  /// \return `true` if default color mode is used, `false` otherwise.
  bool use_default_color() const;

  /// \brief Checks if the default color mode for normals is used.
  /// \return `true` if default color mode for normals is used, `false` otherwise.
  bool use_default_color_normal() const;

  /// \brief Checks if normals are reversed.
  /// \return `true` if normals are reversed, `false` otherwise.
  bool reverse_normal() const;

  /// \brief Checks if the application state is printed within the console.
  /// \return `true` if the application state is printed, `false` otherwise.
  bool print_application_state() const;

  /// \brief Checks if the clipping plane is enabled.
  /// \return `true` if the clipping plane is enabled, `false` otherwise.
  bool clipping_plane_enabled() const;

  /// \brief Checks if the camera is in orthographic mode.
  /// \return `true` if the camera is in orthographic mode, `false` otherwise.
  bool is_orthographic() const;

  /// \brief Checks if the camera is in orthographic mode and the graphics scene is two-dimensional.
  /// \return `true` if the camera is in orthographic mode and the scene is 2D, `false` otherwise.
  bool is_two_dimensional() const;

  /// \brief Gets the clipping plane when enabled.
  /// \return The clipping plane as a `CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3` object.
  CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 clipping_plane() const;

  /// \brief Gets the graphics scene of the viewer.
  /// \return A reference to the `Graphics_scene` object.
  const Graphics_scene& graphics_scene() const;

  /// @}
};

//------------------------------------------------------------------------------
/*!
 \ingroup PkgBasicViewerClasses

This function opens a new window and draws the given `Graphics_scene` (which must have been populated beforehand). The `title` parameter allows specifying a title for the window. This function is blocking, meaning that the program execution will continue as soon as the user closes the window. This function requires `CGAL_GLFW`, and is only available if the macro `CGAL_USE_BASIC_VIEWER_GLFW` is defined. Linking with the CMake target `CGAL::CGAL_Basic_viewer_GLFW` will automatically link with `CGAL_GLFW` and add the definition `CGAL_USE_BASIC_VIEWER_GLFW`.

\param graphic_scene The graphics scene to be drawn.
\param title The title of the window. Defaults to "CGAL Basic Viewer (GLFW)".
*/
void draw_graphics_scene(const Graphics_scene& graphic_scene,
                         const char *title="CGAL Basic Viewer (GLFW)")
{}

} // End namespace GLFW
} // End namespace CGAL
