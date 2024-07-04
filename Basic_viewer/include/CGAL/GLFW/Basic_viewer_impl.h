#ifndef CGAL_BASIC_VIEWER_GLFW_IMPL_H
#define CGAL_BASIC_VIEWER_GLFW_IMPL_H

// #include "Basic_viewer.h"

namespace CGAL 
{
namespace GLFW 
{

  inline 
  void Basic_viewer::error_callback(int error, const char *description)
  {
    std::cerr << "GLFW returned an error:\n\t" << description << "(" << error << ")\n";
  }

  inline
  GLFWwindow* Basic_viewer::create_window(int width, int height, const char* title, bool hidden)
  {
    // Initialise GLFW
    if (!glfwInit())
    {
      std::cerr << "Could not start GLFW\n";
      exit(EXIT_FAILURE);
    }

    // OpenGL 2.1 with compatibilty
    // if (hidden)
    // {
    //   glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    // }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Enable the GLFW runtime error callback function defined previously.
    glfwSetErrorCallback(error_callback);

    // Set additional window options
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, WINDOW_SAMPLES); // MSAA

    // Create window using GLFW
    GLFWwindow* window = glfwCreateWindow(width, height, title, nullptr, nullptr);

    // Ensure the window is set up correctly
    if (!window)
    {
      std::cerr << "Could not open GLFW window\n";
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    // Let the window be the current OpenGL context and initialise glad
    glfwMakeContextCurrent(window);

    // Initialized GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
      std::cerr << "Failed to initialized GLAD!";
      glfwDestroyWindow(window);
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    // Print various OpenGL information to stdout
    std::cout << glGetString(GL_VENDOR) << ": " << glGetString(GL_RENDERER) << '\n';
    std::cout << "GLFW\t " << glfwGetVersionString() << '\n';
    std::cout << "OpenGL\t " << glGetString(GL_VERSION) << '\n';
    std::cout << "GLSL\t " << glGetString(GL_SHADING_LANGUAGE_VERSION) << "\n\n";

    return window;
  }

  inline
  Basic_viewer::Basic_viewer(
    const Graphics_scene* graphicScene,
    const char* title,
    bool drawVertices,
    bool drawEdges,
    bool drawFaces,
    bool drawRays,
    bool drawLines, 
    bool useMonoColor,
    bool inverseNormal,
    bool flatShading
  ) : m_scene(graphicScene),
    m_title(title),
    m_drawVertices(drawVertices),
    m_drawEdges(drawEdges),
    m_drawFaces(drawFaces),
    m_drawRays(drawRays),
    m_drawLines(drawLines),
    m_useMonoColor(useMonoColor),
    m_inverseNormal(inverseNormal),
    m_flatShading(flatShading)
  {
  }

  inline
  void Basic_viewer::initialize(bool screenshotOnly)
  {
    m_window = create_window(m_windowSize.x(), m_windowSize.y(), m_title, screenshotOnly);

    if (!screenshotOnly)
    {
      // Set event callbacks 
      glfwSetWindowUserPointer(m_window, this);
      glfwSetKeyCallback(m_window, key_callback);
      glfwSetCursorPosCallback(m_window, cursor_callback);
      glfwSetMouseButtonCallback(m_window, mouse_btn_callback);
      glfwSetScrollCallback(m_window, scroll_callback);
      glfwSetFramebufferSizeCallback(m_window, window_size_callback);
    }

    GLint openglMajorVersion, openglMinorVersion;
    glGetIntegerv(GL_MAJOR_VERSION, &openglMajorVersion);
    glGetIntegerv(GL_MINOR_VERSION, &openglMinorVersion);

    if (openglMajorVersion > 4 || openglMajorVersion == 4 && openglMinorVersion >= 3)
    {
      m_isOpengl4_3 = true;
    }

    compile_shaders();
    initialize_camera();
    initialize_buffers();
    initialize_and_load_world_axis();
    initialize_and_load_clipping_plane();
    if (!screenshotOnly)
    {
      initialize_keys_actions();
    }
  }

  inline
  void Basic_viewer::show()
  {
    float elapsedTime = 0.0f;
    float lastFrame = 0.0;
    while (!glfwWindowShouldClose(m_window))
    {
      float currentFrame = static_cast<float>(glfwGetTime());
      float deltaTime = currentFrame - lastFrame;
      lastFrame = currentFrame;

      handle_events(deltaTime);
      render_scene(deltaTime);
      glfwSwapBuffers(m_window);

      print_application_state(elapsedTime, deltaTime);
    }

    clear_application();
  }

  void Basic_viewer::clear_application()
  {
    m_plShader.destroy();
    m_faceShader.destroy();
    m_lineShader.destroy();
    m_planeShader.destroy();
    glDeleteBuffers(NB_GL_BUFFERS, m_vbo);
    glDeleteVertexArrays(NB_VAO_BUFFERS, m_vao);
    glfwDestroyWindow(m_window);
    glfwTerminate();
  }

  inline
  void Basic_viewer::make_screenshot(const std::string& filePath)
  {
    render_scene(0);
    glfwSwapBuffers(m_window);
    capture_screenshot(filePath);
    clear_application();
  }

  void generate_grid(Line_renderer& renderer, const vec3f& color, float size, int nbSubdivisions=10)
  {
    for (unsigned int i = 0; i <= nbSubdivisions; ++i)
    {
      float pos = float(size * (2.0 * i / nbSubdivisions - 1.0));
      renderer.add_line(vec3f(pos, -size, 0.f), vec3f(pos, size, 0.f), color);
      renderer.add_line(vec3f(-size, pos, 0.f), vec3f(size, pos, 0.f), color);
    }
  }

  inline
  void Basic_viewer::compile_shaders()
  {
    const char* FACE_VERTEX = m_isOpengl4_3 ? VERTEX_SOURCE_COLOR : VERTEX_SOURCE_COLOR_COMP;
    const char* FACE_FRAGMENT = m_isOpengl4_3 ? FRAGMENT_SOURCE_COLOR : FRAGMENT_SOURCE_COLOR_COMP;
    const char* PL_VERTEX = m_isOpengl4_3 ? VERTEX_SOURCE_P_L : VERTEX_SOURCE_P_L_COMP;
    const char* PL_FRAGMENT = m_isOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    const char* PLANE_VERTEX = VERTEX_SOURCE_CLIPPING_PLANE;
    const char* PLANE_FRAGMENT = FRAGMENT_SOURCE_CLIPPING_PLANE;

    m_faceShader = Shader::load_shader(FACE_VERTEX, FACE_FRAGMENT, "FACE");
    m_plShader = Shader::load_shader(PL_VERTEX, PL_FRAGMENT, "PL");
    m_planeShader = Shader::load_shader(PLANE_VERTEX, PLANE_FRAGMENT, "PLANE");

    // For world axis and grid 
    const char* LINE_VERTEX = VERTEX_SOURCE_LINE;
    const char* LINE_FRAGMENT = FRAGMENT_SOURCE_LINE;
    m_lineShader = Shader::load_shader(LINE_VERTEX, LINE_FRAGMENT, "LINE");
  }

  inline 
  void Basic_viewer::initialize_camera()
  {
    vec3f pmin(
      m_scene->bounding_box().xmin(),
      m_scene->bounding_box().ymin(),
      m_scene->bounding_box().zmin());

    vec3f pmax(
      m_scene->bounding_box().xmax(),
      m_scene->bounding_box().ymax(),
      m_scene->bounding_box().zmax());

    m_camera.lookat(pmin, pmax);  
  }

  inline 
  void Basic_viewer::initialize_and_load_world_axis()
  {
    {
      m_worldAxisRenderer.initialize_buffers();
      m_worldAxisRenderer.set_width(3.f);
      m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
      m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
      m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitZ(), vec3f(0, 0, 1)); // z-axis
      m_worldAxisRenderer.load_buffers();
    }

    float cameraSize = m_camera.get_size() * 0.5;
    {
      m_XYAxisRenderer.initialize_buffers();
      m_XYAxisRenderer.set_width(5.f);
      m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
      m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
      m_XYAxisRenderer.load_buffers();
    }

    {
      m_XYGridRenderer.initialize_buffers();
      m_XYGridRenderer.set_width(2.f);
      m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitX(), vec3f(.8f, .8f, .8f)); // -x-axis
      m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitY(), vec3f(.8f, .8f, .8f)); // -y-axis
      m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitZ(), vec3f(.8f, .8f, .8f)); // -z-axis

      vec3f color(.8f, .8f, .8f);

      generate_grid(m_XYGridRenderer, color, cameraSize);
      
      m_XYGridRenderer.load_buffers();
    }
  }

  inline
  void Basic_viewer::load_buffer(int i, int location, const std::vector<float>& vector, int dataCount)
  {
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[i]);

    glBufferData(GL_ARRAY_BUFFER, vector.size() * sizeof(float), vector.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(location, dataCount, GL_FLOAT, GL_FALSE, dataCount * sizeof(float), nullptr);

    glEnableVertexAttribArray(location);
  }

  inline
  void Basic_viewer::load_buffer(int i, int location, int gsEnum, int dataCount)
  {
    const auto& vector = m_scene->get_array_of_index(gsEnum);
    load_buffer(i, location, vector, dataCount);
  }

  inline
  void Basic_viewer::initialize_buffers()
  {
    glGenBuffers(NB_GL_BUFFERS, m_vbo);
    glGenVertexArrays(NB_VAO_BUFFERS, m_vao);
  }

  inline
  void Basic_viewer::load_scene()
  {
    unsigned int bufn = 0;

    // 1) POINT SHADER

    // 1.1) Mono points
    glBindVertexArray(m_vao[VAO_MONO_POINTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_POINTS, 3);

    // 1.2) Color points
    glBindVertexArray(m_vao[VAO_COLORED_POINTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_POINTS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_POINTS, 3);

    // 2) SEGMENT SHADER

    // 2.1) Mono segments
    glBindVertexArray(m_vao[VAO_MONO_SEGMENTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_SEGMENTS, 3);

    // 2.2) Colored segments
    glBindVertexArray(m_vao[VAO_COLORED_SEGMENTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_SEGMENTS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_SEGMENTS, 3);

    // 3) RAYS SHADER

    // 2.1) Mono segments
    glBindVertexArray(m_vao[VAO_MONO_RAYS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_RAYS, 3);

    // 2.2) Colored segments
    glBindVertexArray(m_vao[VAO_COLORED_RAYS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_RAYS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_RAYS, 3);

    // 4) LINES SHADER

    // 2.1) Mono lines
    glBindVertexArray(m_vao[VAO_MONO_LINES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_LINES, 3);

    // 2.2) Colored lines
    glBindVertexArray(m_vao[VAO_COLORED_LINES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_LINES, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_LINES, 3);

    // 5) FACE SHADER

    // 5.1) Mono faces
    glBindVertexArray(m_vao[VAO_MONO_FACES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_FACES, 3);
    if (m_flatShading)
    {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_MONO_FACES, 3);
    }
    else
    {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_MONO_FACES, 3);
    }

    // 5.2) Colored faces
    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_FACES, 3);
    if (m_flatShading)
    {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_COLORED_FACES, 3);
    }
    else
    {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_COLORED_FACES, 3);
    }
    load_buffer(bufn++, 2, Graphics_scene::COLOR_FACES, 3);

    m_areBuffersInitialized = true;
  }

  inline 
  CGAL::Plane_3<Basic_viewer::Local_kernel> Basic_viewer::clipping_plane() const
  {
    mat4f CPM = m_clippingPlane.get_matrix();
    CGAL::Aff_transformation_3<Basic_viewer::Local_kernel> aff(
      CPM(0, 0), CPM(0, 1), CPM(0, 2), CPM(0, 3),
      CPM(1, 0), CPM(1, 1), CPM(1, 2), CPM(1, 3),
      CPM(2, 0), CPM(2, 1), CPM(2, 2), CPM(2, 3)
    );

    CGAL::Plane_3<Local_kernel> p3(0, 0, 1, 0);
    return p3.transform(aff);
  }

  inline 
  void Basic_viewer::compute_model_view_projection_matrix(const float deltaTime)
  {
    m_camera.update(deltaTime);
    m_clippingPlane.update(deltaTime);

    if (m_animationController.is_running())
    {
      AnimationKeyFrame animationFrame = m_animationController.run();
      m_camera.set_orientation(animationFrame.orientation);
      m_camera.set_position(animationFrame.position);
    }

    m_modelViewMatrix = m_camera.view();

    mat4f projection = m_camera.projection(m_windowSize.x(), m_windowSize.y());
    m_modelViewProjectionMatrix = projection * m_modelViewMatrix;
  }

  inline
  void Basic_viewer::update_uniforms(const float deltaTime)
  {
    compute_model_view_projection_matrix(deltaTime);

    // ================================================================

    update_face_uniforms();
    update_pl_uniforms();
    update_clipping_uniforms();
  }

  inline
  void Basic_viewer::update_face_uniforms()
  {
    m_faceShader.use();

    m_faceShader.set_mat4f("mvp_matrix", m_modelViewProjectionMatrix.data());
    m_faceShader.set_mat4f("mv_matrix", m_modelViewMatrix.data());

    m_faceShader.set_vec4f("light_pos", m_lightPosition.data());
    m_faceShader.set_vec4f("light_diff", m_diffuseColor.data());
    m_faceShader.set_vec4f("light_spec", m_specularColor.data());
    m_faceShader.set_vec4f("light_amb", m_ambientColor.data());
    m_faceShader.set_float("spec_power", m_shininess);

    m_faceShader.set_vec4f("clipPlane", m_clipPlane.data());
    m_faceShader.set_vec4f("pointPlane", m_pointPlane.data());
    m_faceShader.set_float("rendering_transparency", m_clippingPlane.get_transparency());
  }

  inline
  void Basic_viewer::update_pl_uniforms()
  {
    m_plShader.use();

    m_plShader.set_vec4f("clipPlane", m_clipPlane.data());
    m_plShader.set_vec4f("pointPlane", m_pointPlane.data());
    m_plShader.set_mat4f("mvp_matrix", m_modelViewProjectionMatrix.data());
    m_plShader.set_float("point_size", m_sizePoints);
  }

  inline
  void Basic_viewer::update_clipping_uniforms()
  {
    mat4f clippingMatrix = m_clippingPlane.get_matrix();

    m_pointPlane = clippingMatrix * vec4f(0, 0, 0, 1);
    m_clipPlane = clippingMatrix * vec4f(0, 0, 1, 0);

    m_planeShader.use();
    m_planeShader.set_mat4f("vp_matrix", m_modelViewProjectionMatrix.data());
    m_planeShader.set_mat4f("m_matrix", clippingMatrix.data());
  }

  inline
  void Basic_viewer::set_world_axis_uniforms()
  {
    int w = m_windowSize.x();
    int h = m_windowSize.y();

    mat4f view = m_modelViewMatrix;

    // we only want the rotation part of the view matrix  
    mat3f rotation = view.block<3,3>(0,0);

    mat4f rotation4x4 = mat4f::Identity();
    rotation4x4.block<3,3>(0,0) = rotation;

    float aspect = static_cast<float>(w) / h;
    float halfWidth = aspect * 0.1f;
    float halfHeight = 0.1f;
    mat4f proj = ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0f, 1.0f);

    mat4f translate = transform::translation(vec3f(halfWidth - 0.1f*aspect, halfHeight - 0.1f, 0.0f));

    mat4f mvp = proj * rotation4x4 * translate;
    m_lineShader.set_mat4f("mvp_matrix", mvp.data()); 
  }

  inline
  void Basic_viewer::set_XY_grid_uniforms()
  {
    m_lineShader.set_mat4f("mvp_matrix", m_modelViewProjectionMatrix.data());
  }

  inline
  void Basic_viewer::render_scene(const float deltaTime)
  {
    if (!m_areBuffersInitialized)
    {
      load_scene();
    }

    glClearColor(1.0f, 1.0f, 1.0f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);

    // update_uniforms(deltaTime);

    compute_model_view_projection_matrix(deltaTime);

    bool half = m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;

    RenderingMode mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    update_face_uniforms();
    if (m_drawFaces)
    {
      draw_faces();
    }

    update_pl_uniforms();
    if (m_drawVertices)
    {
      draw_vertices(mode);
    }
    if (m_drawEdges)
    {
      draw_edges(mode);
    }
    if (m_drawRays)
    {
      draw_rays();
    }
    if (m_drawLines)
    {
      draw_lines();
    }

    if (m_drawWorldAxis) 
    {
      draw_world_axis();
    }
    if (m_drawXYGrid) 
    {
      draw_xy_grid();
    }
  }

  inline
  vec4f Basic_viewer::color_to_vec4(const CGAL::IO::Color& c) const
  {
    return { static_cast<float>(c.red()) / 255, static_cast<float>(c.green()) / 255, static_cast<float>(c.blue()) / 255, 1.0f };
  }

  inline
  void Basic_viewer::draw_faces()
  {
    if (m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF)
    {
      // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
      // Before rendering all transparent objects, disable z-testing first.

      // 1. draw solid first
      draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

      // 2. draw transparent layer second with back face culling to avoid messy triangles
      glDepthMask(false); // disable z-testing
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      glFrontFace(GL_CW);
      draw_faces_bis(RenderingMode::DRAW_OUTSIDE_ONLY);

      // 3. draw solid again without culling and blend to make sure the solid mesh is visible
      glDepthMask(true); // enable z-testing
      glDisable(GL_CULL_FACE);
      glDisable(GL_BLEND);
      draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

      // 4. render clipping plane here
      render_clipping_plane();
    } 
    else // Not CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF
    {
      if (m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
          m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        // 1. draw solid HALF
        draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

        // 2. render clipping plane here
        render_clipping_plane();
      } 
      else
      {
        // 1. draw solid FOR ALL
        draw_faces_bis(RenderingMode::DRAW_ALL);
      } 
    }
  }

  inline
  void Basic_viewer::draw_faces_bis(RenderingMode mode)
  {
    m_faceShader.set_float("rendering_mode", static_cast<float>(mode));

    vec4f color = color_to_vec4(m_facesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_FACES]);
    glVertexAttrib4fv(2, color.data());
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_FACES));

    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(2);
    }
    else
    {
      glEnableVertexAttribArray(2);
    }
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_FACES));

  }

  inline
  void Basic_viewer::draw_rays()
  {
    m_plShader.set_float("rendering_mode", static_cast<float>(RenderingMode::DRAW_ALL));

    vec4f color = color_to_vec4(m_raysMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_RAYS]);
    glVertexAttrib4fv(1, color.data());

    glLineWidth(m_sizeRays);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_RAYS));

    glBindVertexArray(m_vao[VAO_COLORED_RAYS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_RAYS));
  }

  inline
  void Basic_viewer::draw_vertices(RenderingMode mode)
  {
    m_plShader.set_float("rendering_mode", static_cast<float>(mode));

    vec4f color = color_to_vec4(m_verticeMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_POINTS]);
    glVertexAttrib4fv(1, color.data());
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_POINTS));

    glBindVertexArray(m_vao[VAO_COLORED_POINTS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_POINTS));
  }

  inline
  void Basic_viewer::draw_lines()
  {
    m_plShader.set_float("rendering_mode", static_cast<float>(RenderingMode::DRAW_ALL));

    vec4f color = color_to_vec4(m_linesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_LINES]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_sizeLines);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_LINES));

    glBindVertexArray(m_vao[VAO_COLORED_LINES]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_LINES));
  }

  inline
  void Basic_viewer::draw_edges(RenderingMode mode)
  {
    m_plShader.set_float("rendering_mode", static_cast<float>(mode));

    vec4f color = color_to_vec4(m_edgesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_SEGMENTS]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_sizeEdges);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_SEGMENTS));

    glBindVertexArray(m_vao[VAO_COLORED_SEGMENTS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_SEGMENTS));
  }

  inline
  void Basic_viewer::initialize_and_load_clipping_plane()
  {
    float size = ((m_scene->bounding_box().xmax() - m_scene->bounding_box().xmin()) +
                   (m_scene->bounding_box().ymax() - m_scene->bounding_box().ymin()) +
                   (m_scene->bounding_box().zmax() - m_scene->bounding_box().zmin()));

    const unsigned int NB_SUBDIVISIONS = 30;

    vec3f color(0,0,0);
    m_clippingPlaneRenderer.initialize_buffers();
    m_clippingPlaneRenderer.set_width(0.1f);
    generate_grid(m_clippingPlaneRenderer, color, size, NB_SUBDIVISIONS);
    m_clippingPlaneRenderer.load_buffers();

    m_clippingPlane.set_size(m_camera.get_size());  
  }

  inline
  void Basic_viewer::render_clipping_plane()
  {
    if (!m_drawClippingPlane || !m_isOpengl4_3)
    {
      return;
    }

    update_clipping_uniforms();
    m_clippingPlaneRenderer.draw();
  }

  inline
  void Basic_viewer::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    Basic_viewer* viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_key_event(key, scancode, action, mods);
  }

  inline
  void Basic_viewer::cursor_callback(GLFWwindow* window, double xpos, double ypo)
  {
    Basic_viewer* viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_cursor_event(xpos, ypo);
  }

  inline
  void Basic_viewer::mouse_btn_callback(GLFWwindow* window, int button, int action, int mods)
  {
    Basic_viewer* viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_mouse_button_event(button, action, mods);
  }

  inline
  void Basic_viewer::window_size_callback(GLFWwindow* window, int width, int height)
  {
    Basic_viewer* viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));

    viewer->m_windowSize = {width, height};

    glViewport(0, 0, width, height);
  }

  inline
  void Basic_viewer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
  {
    Basic_viewer* viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_scroll_event(xoffset, yoffset);
  }

  inline 
  void Basic_viewer::print_application_state(float& elapsedTime, const float deltaTime)
  {
    elapsedTime += deltaTime;
    if (elapsedTime * 1000 > 100) // update terminal display each 100ms
    {
      elapsedTime = 0.0f;
      if (m_printApplicationState)
      {
        std::cout << "\33[2K"  
                  << "FPS: "                        << std::round(1 / deltaTime)                        << "\n\33[2K" 
                  << "Camera translation speed: "   << m_camera.get_translation_speed()                 << "    " 
                  << "Camera rotation speed: "      << std::round(m_camera.get_rotation_speed())        << "    "
                  << "CP translation speed: "       << m_clippingPlane.get_translation_speed()          << "\n\33[2K"            
                  << "CP rotation speed: "          << std::round(m_clippingPlane.get_rotation_speed()) << "    "     
                  << "CP constraint axis: "         << m_clippingPlane.get_constraint_axis()            << "\n\33[2K"     
                  << "Light color: ("               << m_ambientColor.x()                               << ", " 
                                                    << m_ambientColor.y()                               << ", " 
                                                    << m_ambientColor.z()                               << ")   "     
                  << "\033[F\033[F\033[F\r" << std::flush;
      }
    }
  }

  inline
  void Basic_viewer::start_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case ROTATE_CLIPPING_PLANE:
    case TRANSLATE_CLIPPING_PLANE:
    case TRANSLATE_CP_IN_CAMERA_DIRECTION:
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      // glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      break;
    }
  }

  inline 
  void Basic_viewer::exit_app() 
  {
    clear_application();

    std::cout << "\n APPLICATION EXITED" << std::endl;
    exit(EXIT_SUCCESS);
  }

    inline
  void Basic_viewer::double_click_event(int btn)
  {
    if (m_camera.is_orbiter()) 
    {
      if (btn == GLFW_MOUSE_BUTTON_RIGHT)
      {
        if (is_key_pressed(m_window, GLFW_KEY_LEFT_CONTROL)) 
        {
          m_camera.reset_orientation();
        }
        m_camera.reset_position();
      }
      if (btn == GLFW_MOUSE_BUTTON_LEFT)
      {
        m_camera.align_to_nearest_axis();
      }
    }
  }

  inline
  void Basic_viewer::scroll_event(const float deltaTime)
  {
    float yoffset = get_scroll_yOffset();

    if (is_key_pressed(m_window, GLFW_KEY_LEFT_SHIFT))
    {
      m_camera.increase_zoom_smoothness(yoffset);
    }
    else if (is_key_pressed(m_window, GLFW_KEY_Z) && !m_camera.is_orthographic())
    {
      m_camera.increase_fov(yoffset);
      m_clippingPlane.set_size(m_camera.get_size());
    }
    else if (is_key_pressed(m_window, GLFW_KEY_LEFT_CONTROL) && m_displayMode != DisplayMode::CLIPPING_PLANE_OFF)
    {
        m_clippingPlane.translation(8.f * yoffset * deltaTime);
    }
    else 
    {
      m_camera.move(8.f * yoffset * deltaTime);
      m_clippingPlane.set_size(m_camera.get_size());
    }
  }

  inline 
  void Basic_viewer::action_event(int action, const float deltaTime)
  {
    switch (action)
    {
    /*APPLICATION*/
    case EXIT: 
      exit_app();
    case PRINT_APPLICATION_STATE:
      m_printApplicationState = !m_printApplicationState; 
      break;
    /*WINDOW*/
    case FULLSCREEN:
      fullscreen();
      break;
    case SCREENSHOT:
      capture_screenshot("./screenshot.png");
      break;
    /*SCENE*/
    case VERTICES_DISPLAY:
      m_drawVertices = !m_drawVertices;
      break;
    case FACES_DISPLAY:
      m_drawFaces = !m_drawFaces;
      break;
    case EDGES_DISPLAY:
      m_drawEdges = !m_drawEdges;
      break;
    case SHADING_MODE:
      m_flatShading = !m_flatShading;
      m_areBuffersInitialized = false;
      break;
    case INVERSE_NORMAL:
      m_inverseNormal = !m_inverseNormal;
      m_scene->reverse_all_normals();
      m_areBuffersInitialized = false;
      break;
    case MONO_COLOR:
      m_useMonoColor = !m_useMonoColor;
      break;
    case INC_EDGES_SIZE:
      m_sizeEdges = std::min(50.f, m_sizeEdges + 3.0f*deltaTime);
      break;
    case DEC_EDGES_SIZE:
      m_sizeEdges = std::max(1.f, m_sizeEdges - 3.0f*deltaTime);
      break;
    case INC_POINTS_SIZE:
      m_sizePoints = std::min(30.f, m_sizePoints + 3.0f*deltaTime);
      break;
    case DEC_POINTS_SIZE:
      m_sizePoints = std::max(1.f, m_sizePoints - 3.0f*deltaTime);
      break;
    case INC_LIGHT_ALL:
      increase_light_all(deltaTime);
      break;
    case DEC_LIGHT_ALL:
      increase_light_all(-deltaTime);
      break;
    case INC_LIGHT_R:
      increase_red_component(deltaTime);
      break;
    case INC_LIGHT_G:
      increase_green_component(deltaTime);
      break;
    case INC_LIGHT_B:
      increase_blue_component(deltaTime);
      break;
    case DEC_LIGHT_R:
      increase_red_component(-deltaTime);
      break;
    case DEC_LIGHT_G:
      increase_green_component(-deltaTime);
      break;
    case DEC_LIGHT_B:
      increase_blue_component(-deltaTime);
      break;
    case DISPLAY_WORLD_AXIS:
      m_drawWorldAxis = !m_drawWorldAxis;
      break;
    case DISPLAY_XY_GRID:
      m_drawXYGrid = !m_drawXYGrid;
      break;
    /*CAMERA*/
    case UP:
      m_camera.move_up(deltaTime);
      break;
    case DOWN:
      m_camera.move_down(deltaTime);
      break;
    case LEFT:
      m_camera.move_left(deltaTime);
      break;
    case RIGHT:
      m_camera.move_right(deltaTime);
      break;
    case FORWARD:
      m_camera.move(deltaTime);
      break;
    case BACKWARDS:
      m_camera.move(-deltaTime);
      break;
    case INC_CAMERA_ROTATION_SMOOTHNESS:
      m_camera.increase_rotation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SMOOTHNESS:
      m_camera.decrease_rotation_smoothness(deltaTime);
      break;
    case INC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_camera.increase_translation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_camera.decrease_translation_smoothness(deltaTime);
      break;
    case SWITCH_CAM_MODE:
      m_camera.toggle_mode();
      break;
    case SWITCH_CAM_TYPE:
      m_camera.toggle_type();
      break;
    case INC_CAMERA_TRANSLATION_SPEED:
      m_camera.increase_translation_speed(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SPEED:
      m_camera.decrease_translation_speed(deltaTime);
      break;
    case INC_CAMERA_ROTATION_SPEED:
      m_camera.increase_rotation_speed(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SPEED:
      m_camera.decrease_rotation_speed(deltaTime);
      break;
    case ROTATE_CAMERA:
      rotate_camera();
      break;
    case TRANSLATE_CAMERA:
      translate_camera(deltaTime);
      break;
    case RESET_CAM:
      reset_camera_and_clipping_plane();
      break;
    /*CLIPPING PLANE*/
    case INC_CP_TRANSLATION_SPEED:
      m_clippingPlane.increase_translation_speed(deltaTime);
      break;
    case DEC_CP_TRANSLATION_SPEED:
      m_clippingPlane.decrease_translation_speed(deltaTime);
      break;
    case INC_CP_ROTATION_SPEED:
      m_clippingPlane.increase_rotation_speed(deltaTime);
      break;
    case DEC_CP_ROTATION_SPEED:
      m_clippingPlane.decrease_rotation_speed(deltaTime);
      break;
    case TOGGLE_CLIPPING_PLANE_DISPLAY:
      m_drawClippingPlane = !m_drawClippingPlane;
      break;
    case TOGGLE_CLIPPING_PLANE_MODE:
      switch_display_mode();
      break;
    case ROTATE_CLIPPING_PLANE:
      rotate_clipping_plane();
      break;
    case TRANSLATE_CLIPPING_PLANE:
      translate_clipping_plane(deltaTime);
      break;
    case TRANSLATE_CP_IN_CAMERA_DIRECTION:
      translate_clipping_plane(deltaTime, true);
      break;
    case TOGGLE_CP_CONSTRAINT_AXIS:
      m_clippingPlane.switch_constraint_axis();
      break;
    /*ANIMATION*/
    case SAVE_KEY_FRAME:
      save_key_frame();
      break;
    case RUN_OR_STOP_ANIMATION:
      run_or_stop_animation();
      break;
    case CLEAR_ANIMATION: 
      m_animationController.clear_buffer();
      break;
    }
  }

  inline 
  void Basic_viewer::increase_red_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.x() += speed;
    if (m_ambientColor.x() > 1.f)
      m_ambientColor.x() = 1.f;
    if (m_ambientColor.x() < 0.f)
      m_ambientColor.x() = 0.f;

    m_diffuseColor.x() += speed;
    if (m_diffuseColor.x() > 1.f)
      m_diffuseColor.x() = 1.f;
    if (m_diffuseColor.x() < 0.f)
      m_diffuseColor.x() = 0.f;

    m_specularColor.x() += speed;
    if (m_specularColor.x() > 1.f)
      m_specularColor.x() = 1.f;
    if (m_specularColor.x() < 0.f)
      m_specularColor.x() = 0.f;
  }

  inline 
  void Basic_viewer::increase_green_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.y() += speed;
    if (m_ambientColor.y() > 1.f)
      m_ambientColor.y() = 1.f;
    if (m_ambientColor.y() < 0.f)
      m_ambientColor.y() = 0.f;

    m_diffuseColor.y() += speed;
    if (m_diffuseColor.y() > 1.f)
      m_diffuseColor.y() = 1.f;
    if (m_diffuseColor.y() < 0.f)
      m_diffuseColor.y() = 0.f;

    m_specularColor.y() += speed;
    if (m_specularColor.y() > 1.f)
      m_specularColor.y() = 1.f;
    if (m_specularColor.y() < 0.f)
      m_specularColor.y() = 0.f;
  }

  inline 
  void Basic_viewer::increase_blue_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.z() += speed;
    if (m_ambientColor.z() > 1.f)
      m_ambientColor.z() = 1.f;
    if (m_ambientColor.z() < 0.f)
      m_ambientColor.z() = 0.f;

    m_diffuseColor.z() += speed;
    if (m_diffuseColor.z() > 1.f)
      m_diffuseColor.z() = 1.f;
    if (m_diffuseColor.z() < 0.f)
      m_diffuseColor.z() = 0.f;

    m_specularColor.z() += speed;
    if (m_specularColor.z() > 1.f)
      m_specularColor.z() = 1.f;
    if (m_specularColor.z() < 0.f)
      m_specularColor.z() = 0.f;
  }


  inline 
  void Basic_viewer::increase_light_all(const float deltaTime)
  {
    increase_red_component(deltaTime);
    increase_green_component(deltaTime);
    increase_blue_component(deltaTime);
  }

  inline 
  void Basic_viewer::save_key_frame() 
  {
    if (m_camera.is_orbiter())
    {
      m_animationController.add_key_frame(
        m_camera.get_position(),
        m_camera.get_orientation()
      );
    }
  }

  inline 
  void Basic_viewer::run_or_stop_animation() 
  {
    if (m_camera.is_orbiter())
    {
      if (m_animationController.is_running()) 
      {
        m_animationController.stop(m_animationController.get_frame());
      }
      else 
      {
        m_animationController.start();
      }
    }
  }

  inline
  void Basic_viewer::reset_camera_and_clipping_plane()
  {
    m_camera.reset_all();
    m_clippingPlane.reset_all();
    m_clippingPlane.set_size(m_camera.get_size());
  }

  inline
  void Basic_viewer::switch_display_mode()
  {
    switch(m_displayMode) 
    {
    case DisplayMode::CLIPPING_PLANE_OFF:
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF:
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF: 
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY: 
      m_displayMode = DisplayMode::CLIPPING_PLANE_OFF; 
      break;
    }
  }

  inline
  void Basic_viewer::rotate_clipping_plane()
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();
    
    m_clippingPlane.rotation(
      mouseDeltaX, 
      mouseDeltaY
    );
  }

  inline
  void Basic_viewer::translate_clipping_plane(const float deltaTime, bool useCameraForward)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    float deltaX = deltaTime * mouseDeltaX;
    float deltaY = deltaTime * mouseDeltaY;

    if (useCameraForward)
    {
      vec3f forwardDirection = m_camera.get_forward();

    float s = deltaX;
    if (abs(deltaY) > abs(deltaX))
      s = -deltaY;

      m_clippingPlane.translation(forwardDirection, s);
    }
    else 
    {
      m_clippingPlane.translation(-deltaX, deltaY);
    }
  }

  inline
  void Basic_viewer::rotate_camera()
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();
    
    m_camera.rotation(
      mouseDeltaX, 
      mouseDeltaY
    );

    m_clippingPlane.set_up_axis(m_camera.get_up());
    m_clippingPlane.set_right_axis(m_camera.get_right());
  }

  inline
  void Basic_viewer::translate_camera(const float deltaTime)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    m_camera.translation(
      deltaTime * -mouseDeltaX,
      deltaTime * mouseDeltaY
    );
  }

  inline
  void Basic_viewer::fullscreen()
  {
    m_isFullscreen = !m_isFullscreen;

    if (m_isFullscreen)
    {
      int count;
      m_oldWindowSize = m_windowSize;

      GLFWmonitor *monitor = glfwGetMonitors(&count)[0];
      const GLFWvidmode *mode = glfwGetVideoMode(monitor);

      glfwGetWindowPos(m_window, &m_oldWindowPosition.x(), &m_oldWindowPosition.y());
      glfwSetWindowMonitor(m_window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
      glViewport(0, 0, mode->width, mode->height);

      m_windowSize.x() = mode->width;
      m_windowSize.y() = mode->height;
      return;
    }

    m_windowSize = m_oldWindowSize;
    glfwSetWindowMonitor(m_window, nullptr, m_oldWindowPosition.x(), m_oldWindowPosition.y(), m_windowSize.x(), m_windowSize.y(), 60);
    glViewport(0, 0, m_windowSize.x(), m_windowSize.y());
  }

  inline
  void Basic_viewer::capture_screenshot(const std::string &filepath)
  {
    // https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/ (thanks)
    // https://github.com/nothings/stb/
    // The stb lib used here is from glfw/deps

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);

    const GLsizei NB_CHANNELS = 4;
    GLsizei stride = NB_CHANNELS * m_windowSize.x();
    stride += (stride % 4) ? (4 - stride % 4) : 0; // stride must be a multiple of 4
    GLsizei bufferSize = stride * m_windowSize.y();

    std::vector<char> buffer(bufferSize);
    m_faceShader.use(); 
    glPixelStorei(GL_PACK_ALIGNMENT, 4);

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, m_windowSize.x(), m_windowSize.y(), GL_RGBA, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath.data(), m_windowSize.x(), m_windowSize.y(), NB_CHANNELS, buffer.data(), stride);
  }

  void Basic_viewer::draw_world_axis() 
  {
    int &w = m_windowSize.x();
    int &h = m_windowSize.y();

    set_world_axis_uniforms();

    glViewport(w - w / 5, h - h / 5, w / 5, h / 5);
    m_worldAxisRenderer.draw();

    // Restore the main viewport
    glViewport(0, 0, w, h);
  }

  void Basic_viewer::draw_xy_grid() 
  { 
    set_XY_grid_uniforms();
    m_XYAxisRenderer.draw();
    m_XYGridRenderer.draw();
  }
 
  inline
  void Basic_viewer::end_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
      break;
    }
  }

  inline
  void Basic_viewer::initialize_keys_actions()
  {
    set_keyboard_layout(KeyboardLayout::AZERTY);

    /*APPLICATION*/
    add_keyboard_action({GLFW_KEY_ESCAPE},                   InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_Q, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_T},                        InputMode::RELEASE, PRINT_APPLICATION_STATE);

    /*WINDOW*/
    add_keyboard_action({GLFW_KEY_ENTER, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, FULLSCREEN);
    add_keyboard_action({GLFW_KEY_F2},                       InputMode::RELEASE, SCREENSHOT);

    /*SCENE*/
    add_keyboard_action({GLFW_KEY_A}, InputMode::RELEASE, DISPLAY_WORLD_AXIS);
    add_keyboard_action({GLFW_KEY_G}, InputMode::RELEASE, DISPLAY_XY_GRID);

    add_keyboard_action({GLFW_KEY_W}, InputMode::RELEASE, FACES_DISPLAY);
    add_keyboard_action({GLFW_KEY_V}, InputMode::RELEASE, VERTICES_DISPLAY);
    add_keyboard_action({GLFW_KEY_E}, InputMode::RELEASE, EDGES_DISPLAY);

    add_keyboard_action({GLFW_KEY_S}, InputMode::RELEASE, SHADING_MODE);
    add_keyboard_action({GLFW_KEY_N}, InputMode::RELEASE, INVERSE_NORMAL);
    add_keyboard_action({GLFW_KEY_M}, InputMode::RELEASE, MONO_COLOR);

    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD, GLFW_KEY_LEFT_CONTROL},                     InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_6, GLFW_KEY_LEFT_CONTROL},                          InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT, GLFW_KEY_LEFT_CONTROL},                InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT},                        InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD},                                            InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_6},                                                 InputMode::HOLD, DEC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT},                                       InputMode::HOLD, DEC_EDGES_SIZE);

    add_keyboard_action({GLFW_KEY_PAGE_UP},                              InputMode::HOLD, INC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN},                            InputMode::HOLD, DEC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_SHIFT},      InputMode::HOLD, INC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_SHIFT},      InputMode::HOLD, DEC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_ALT},        InputMode::HOLD, INC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_ALT},        InputMode::HOLD, DEC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_CONTROL},    InputMode::HOLD, INC_LIGHT_B);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_CONTROL},    InputMode::HOLD, DEC_LIGHT_B);

    /*CAMERA*/
    add_keyboard_action({GLFW_KEY_UP,   GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, FORWARD);
    add_keyboard_action({GLFW_KEY_DOWN, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, BACKWARDS);

    add_keyboard_action({GLFW_KEY_UP},    InputMode::HOLD, UP);
    add_keyboard_action({GLFW_KEY_DOWN},  InputMode::HOLD, DOWN);
    add_keyboard_action({GLFW_KEY_LEFT},  InputMode::HOLD, LEFT);
    add_keyboard_action({GLFW_KEY_RIGHT}, InputMode::HOLD, RIGHT);

    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, RESET_CAM);
    
    add_keyboard_action({GLFW_KEY_SPACE}, InputMode::RELEASE, SWITCH_CAM_TYPE);
    add_keyboard_action({GLFW_KEY_O},     InputMode::RELEASE, SWITCH_CAM_MODE); 

    add_keyboard_action({GLFW_KEY_X},                         InputMode::HOLD, INC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R},                         InputMode::HOLD, INC_CAMERA_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT},    InputMode::HOLD, DEC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT},    InputMode::HOLD, DEC_CAMERA_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT},  InputMode::HOLD, ROTATE_CAMERA);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT}, InputMode::HOLD, TRANSLATE_CAMERA);

    /*CLIPPING PLANE*/
    add_keyboard_action({GLFW_KEY_C, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, TOGGLE_CLIPPING_PLANE_DISPLAY);
    add_keyboard_action({GLFW_KEY_C},                    InputMode::RELEASE, TOGGLE_CLIPPING_PLANE_MODE);

    add_keyboard_action({GLFW_KEY_A, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, TOGGLE_CP_CONSTRAINT_AXIS);

    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL},                      InputMode::HOLD, INC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL},                      InputMode::HOLD, INC_CP_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT,   GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, ROTATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT,  GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_MIDDLE, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CP_IN_CAMERA_DIRECTION);

    /*ANIMATION*/
    add_keyboard_action({GLFW_KEY_F1},                                InputMode::RELEASE, RUN_OR_STOP_ANIMATION);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_LEFT_ALT},             InputMode::RELEASE, SAVE_KEY_FRAME);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_D, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, CLEAR_ANIMATION);

    /*===================== BIND DESCRIPTIONS ============================*/

    set_action_description({
      {"Animation", {
        {get_binding_text_from_action(SAVE_KEY_FRAME),        "Add a key frame to animation (only available in orbiter)"},
        {get_binding_text_from_action(CLEAR_ANIMATION),       "Delete animation (only available in orbiter)"},
        {get_binding_text_from_action(RUN_OR_STOP_ANIMATION), "Start/Stop animation (only available in orbiter)"},
      }},
      {"Clipping plane", {
        {get_binding_text_from_action(TOGGLE_CLIPPING_PLANE_MODE),        "Switch clipping plane display mode"},
        {get_binding_text_from_action(TOGGLE_CLIPPING_PLANE_DISPLAY),     "Toggle clipping plane rendering on/off"},
        {get_binding_text_from_action(TOGGLE_CP_CONSTRAINT_AXIS),         "Toggle constraint axis for clipping plane rotation"},
        {"",                                                              ""},
        {get_binding_text_from_action(INC_CP_ROTATION_SPEED),             "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CP_ROTATION_SPEED),             "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CP_TRANSLATION_SPEED),          "Increase translation speed"},
        {get_binding_text_from_action(DEC_CP_TRANSLATION_SPEED),          "Decrease translation speed"},
        {"",                                                              ""},
        {get_binding_text_from_action(ROTATE_CLIPPING_PLANE),             "Rotate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CLIPPING_PLANE),          "Translate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CP_IN_CAMERA_DIRECTION),  "Translate the clipping plane along camera direction axis when enabled"},
      }},
      {"Camera", {
        {get_binding_text_from_action(FORWARD),                               "Move camera forward"},
        {get_binding_text_from_action(BACKWARDS),                             "Move camera backwards"},
        {get_binding_text_from_action(UP),                                    "Move camera up"},
        {get_binding_text_from_action(DOWN),                                  "Move camera down"},
        {get_binding_text_from_action(RIGHT),                                 "Move camera right"},
        {get_binding_text_from_action(LEFT),                                  "Move camera left"},
        {"",                                                                  ""},
        {get_binding_text_from_action(SWITCH_CAM_MODE),                       "Switch to Perspective/Orthographic view"},
        {get_binding_text_from_action(SWITCH_CAM_TYPE),                       "Switch to Orbiter/Free-fly camera type"},
        {"",                                                                  ""},
        {get_binding_text_from_action(ROTATE_CAMERA),                         "Rotate the camera"},
        {get_binding_text_from_action(TRANSLATE_CAMERA),                      "Translate the camera"},
        {get_binding_text_from_action(RESET_CAM),                             "Reset camera"},
        {"",                                                                  ""},
        {get_binding_text_from_action(INC_CAMERA_ROTATION_SPEED),             "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CAMERA_ROTATION_SPEED),             "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CAMERA_TRANSLATION_SPEED),          "Increase translation speed"},
        {get_binding_text_from_action(DEC_CAMERA_TRANSLATION_SPEED),          "Decrease translation speed"},
        {"",                                                                  ""},
        {"WHEEL",                                                             "Zoom in/out"},
        {"LEFT_DOUBLE_CLICK",                                                 "Center the camera to the object"},
        {"LCTRL+LEFT_DOUBLE_CLICK",                                           "Reset camera position & orientation"},
        {"RIGHT_DOUBLE_CLICK",                                                "Aligns camera to nearest axis"},
        {"Z+WHEEL",                                                           "Increase/Decrease FOV"},
      }},
      {"Scene", {
        {get_binding_text_from_action(INC_LIGHT_ALL),       "Increase light (all colors, use shift/alt/ctrl for one rgb component)"},
        {get_binding_text_from_action(DEC_LIGHT_ALL),       "Decrease light (all colors, use shift/alt/ctrl for one rgb component)"},
        {"",                                                ""},
        {get_binding_text_from_action(VERTICES_DISPLAY),    "Toggle vertices display"},
        {get_binding_text_from_action(EDGES_DISPLAY),       "Toggle edges display"},
        {get_binding_text_from_action(FACES_DISPLAY),       "Toggle faces display"},
        {"",                                                ""},
        {get_binding_text_from_action(DISPLAY_WORLD_AXIS),  "Toggle the display of the world axis"},
        {get_binding_text_from_action(DISPLAY_XY_GRID),     "Toggle the display of the XY grid"},
        {"",                                                ""},
        {get_binding_text_from_action(INC_POINTS_SIZE),     "Increase size of vertices"},
        {get_binding_text_from_action(DEC_POINTS_SIZE),     "Decrease size of vertices"},
        {get_binding_text_from_action(INC_EDGES_SIZE),      "Increase size of edges"},
        {get_binding_text_from_action(DEC_EDGES_SIZE),      "Decrease size of edges"},
        {"",                                                ""},
        {get_binding_text_from_action(MONO_COLOR),          "Toggle mono color"},
        {get_binding_text_from_action(INVERSE_NORMAL),      "Invert direction of normals"},
        {get_binding_text_from_action(SHADING_MODE),        "Switch between flat/Gouraud shading display"},
      }},         
      {"Window", {
        {get_binding_text_from_action(FULLSCREEN),  "Switch to windowed/fullscreen mode"},
        {get_binding_text_from_action(SCREENSHOT),  "Take a screenshot of the current view"},
      }},
      {"Application", {
        {get_binding_text_from_action(EXIT),                    "Exit program"},
        {"",                                                    ""},
        {get_binding_text_from_action(PRINT_APPLICATION_STATE), "Print application state within the terminal (FPS...)"},
      }},
    });

    print_help();
  }
} // end namespace GLFW 
} // end namespace CGAL

#endif // CGAL_BASIC_VIEWER_GLFW_IMPL_H