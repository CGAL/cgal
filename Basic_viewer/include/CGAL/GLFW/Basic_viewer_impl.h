#pragma once

#include "Basic_viewer.h"

namespace CGAL {
namespace GLFW {
  void Basic_viewer::error_callback(int error, const char *description)
  {
    fprintf(stderr, "GLFW returned an error:\n\t%s (%i)\n", description, error);
  }

  GLFWwindow *Basic_viewer::create_window(int width, int height, const char *title, bool hidden)
  {
    // Initialise GLFW
    if (!glfwInit())
    {
      fprintf(stderr, "Could not start GLFW\n");
      exit(EXIT_FAILURE);
    }

    // OpenGL 2.1 with compatibilty
    if (hidden)
      glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Enable the GLFW runtime error callback function defined previously.
    glfwSetErrorCallback(error_callback);

    // Set additional window options
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, windowSamples); // MSAA

    // Create window using GLFW
    GLFWwindow *window = glfwCreateWindow(width,
                                          height,
                                          title,
                                          nullptr,
                                          nullptr);

    // Ensure the window is set up correctly
    if (!window)
    {
      fprintf(stderr, "Could not open GLFW window\n");
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    // Let the window be the current OpenGL context and initialise glad
    glfwMakeContextCurrent(window);

    // Initialized GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
      std::cerr << "Failed to initialized GLAD!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // gladLoadGL();

    // White background

    // Print various OpenGL information to stdout
    printf("%s: %s\n", glGetString(GL_VENDOR), glGetString(GL_RENDERER));
    printf("GLFW\t %s\n", glfwGetVersionString());
    printf("OpenGL\t %s\n", glGetString(GL_VERSION));
    printf("GLSL\t %s\n\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

    return window;
  }

  Basic_viewer::Basic_viewer(
      const Graphics_scene *graphicScene,
      const char *title,
      bool draw_vertices,
      bool draw_edges,
      bool draw_faces,
      bool use_mono_color,
      bool inverse_normal,
      bool draw_rays,
      bool draw_lines) : m_scene(graphicScene),
                         m_title(title),
                         m_drawVertices(draw_vertices),
                         m_drawEdges(draw_edges),
                         m_drawRays(draw_rays),
                         m_drawLines(draw_lines),
                         m_drawFaces(draw_faces),
                         m_useMonoColor(use_mono_color),
                         m_inverseNormal(inverse_normal)
  {
    init_keys_actions();
  }

  void Basic_viewer::show()
  {
    m_window = create_window(m_window_size.x(), m_window_size.y(), m_title);
    init_buffers();

    glGenBuffers(NB_GL_BUFFERS, m_vbo);
    glGenVertexArrays(NB_VAO_BUFFERS, m_vao);

    glfwSetWindowUserPointer(m_window, this);
    glfwSetKeyCallback(m_window, key_callback);
    glfwSetCursorPosCallback(m_window, cursor_callback);
    glfwSetMouseButtonCallback(m_window, mouse_btn_callback);
    glfwSetScrollCallback(m_window, scroll_callback);
    glfwSetFramebufferSizeCallback(m_window, window_size_callback);

    print_help();
    set_cam_mode(m_camMode);

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major > 4 || major == 4 && minor >= 3)
    {
      m_is_opengl_4_3 = true;
    }

    compile_shaders();

    vec3f pmin(
        m_scene->bounding_box().xmin(),
        m_scene->bounding_box().ymin(),
        m_scene->bounding_box().zmin());

    vec3f pmax(
        m_scene->bounding_box().xmax(),
        m_scene->bounding_box().ymax(),
        m_scene->bounding_box().zmax());

    m_camera.lookat(pmin, pmax);

    double lastFrame = glfwGetTime();
    while (!glfwWindowShouldClose(m_window))
    {
      double currentFrame = glfwGetTime();
      double deltaTime = currentFrame - lastFrame;
      lastFrame = currentFrame;
      handle_events(deltaTime);
      render_scene();
      glfwSwapBuffers(m_window);
      // update_frame();
    }

    glfwTerminate();
  }

  void Basic_viewer::make_screenshot(const std::string &pngpath)
  {
    m_areBuffersInitialized = false;
    m_window = create_window(m_window_size.x(), m_window_size.y(), m_title, true);
    init_buffers();

    set_cam_mode(m_camMode);

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major > 4 || major == 4 && minor >= 3)
    {
      m_is_opengl_4_3 = true;
    }

    compile_shaders();
    render_scene();
    glfwSwapBuffers(m_window);
    screenshot(pngpath);
    glfwTerminate();
  }

  void Basic_viewer::update_frame()
  {
    vec3f dO, dx, dy;
    float z = 0.0f; // Définissez la profondeur souhaitée

    // Obtenir le cadre de référence à la profondeur z
    m_camera.frame(z, dO, dx, dy);

    // Afficher les vecteurs de direction pour le débogage
    std::cout << "\x1B[2J\x1B[H";
    // std::cout << "Origin: " << dO.x() << ", " << dO.y() << ", " << dO.z() << std::endl;
    std::cout << "Forward: " << m_camera.forward().x() << ", " << m_camera.forward().y() << ", " << m_camera.forward().z() << std::endl;
    std::cout << "Direction X: " << dx.x() << ", " << dx.y() << ", " << dx.z() << std::endl;
    std::cout << "Direction Y: " << dy.x() << ", " << dy.y() << ", " << dy.z() << std::endl;
    vec3f position = m_camera.position();
    // vec3f test = multVecMat(m_camera.position(), m_camera.viewport());
    // std::cout << "Viewport : " << test.x() << ", " << test.y() << ", " << test.z() << std::endl;
    std::cout << "Viewport : " << position.x() << ", " << position.y() << ", " << position.z() << std::endl;
  }

  void Basic_viewer::compile_shaders()
  {
    const char *face_vert = m_is_opengl_4_3 ? vertex_source_color : vertex_source_color_comp;
    const char *face_frag = m_is_opengl_4_3 ? fragment_source_color : fragment_source_color_comp;
    const char *pl_vert = m_is_opengl_4_3 ? vertex_source_p_l : vertex_source_p_l_comp;
    const char *pl_frag = m_is_opengl_4_3 ? fragment_source_p_l : fragment_source_p_l_comp;
    const char *plane_vert = vertex_source_clipping_plane;
    const char *plane_frag = fragment_source_clipping_plane;

    m_face_shader = Shader::loadShader(face_vert, face_frag, "FACE");
    m_pl_shader = Shader::loadShader(pl_vert, pl_frag, "PL");
    m_plane_shader = Shader::loadShader(plane_vert, plane_frag, "PLANE");
  }

  void Basic_viewer::load_buffer(int i, int location, const std::vector<float> &vector, int dataCount)
  {
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[i]);

    glBufferData(GL_ARRAY_BUFFER, vector.size() * sizeof(float), vector.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(location, dataCount, GL_FLOAT, GL_FALSE, dataCount * sizeof(float), nullptr);

    glEnableVertexAttribArray(location);
  }

  void Basic_viewer::load_buffer(int i, int location, int gsEnum, int dataCount)
  {
    const std::vector<float> &vector = m_scene->get_array_of_index(gsEnum);
    load_buffer(i, location, vector, dataCount);
  }

  void Basic_viewer::init_buffers()
  {
    if (m_areBuffersInitialized)
    {
      glGenBuffers(NB_GL_BUFFERS, m_vbo);
      glGenVertexArrays(NB_VAO_BUFFERS, m_vao);
      m_areBuffersInitialized = true;
    }
  }

  void Basic_viewer::load_scene()
  {
    unsigned int bufn = 0;

    // 1) POINT SHADER

    // 1.1) Mono points
    m_pl_shader.use();
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
    m_face_shader.use();
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

    // 6) clipping plane shader
    if (m_is_opengl_4_3)
    {
      generate_clipping_plane();
      glBindVertexArray(m_vao[VAO_CLIPPING_PLANE]);
      load_buffer(bufn++, 0, m_array_for_clipping_plane, 3);
    }

    m_areBuffersInitialized = true;
  }

  CGAL::Plane_3<Basic_viewer::Local_kernel> Basic_viewer::clipping_plane() const
  {
    const mat4f cpm = m_clipping_matrix;
    CGAL::Aff_transformation_3<Basic_viewer::Local_kernel> aff(
        cpm(0, 0), cpm(0, 1), cpm(0, 2), cpm(0, 3),
        cpm(1, 0), cpm(1, 1), cpm(1, 2), cpm(1, 3),
        cpm(2, 0), cpm(2, 1), cpm(2, 2), cpm(2, 3));

    CGAL::Plane_3<Local_kernel> p3(0, 0, 1, 0);
    return p3.transform(aff);
  }

  void Basic_viewer::update_uniforms()
  {
    // m_mv = lookAt(m_camPosition, m_camPosition + m_camForward, vec3f(0,1,0)) * m_scene_rotation;

    mat4f view = m_camera.view();
    mat4f projection = m_camera.projection(m_window_size.x(), m_window_size.y());

    m_mv = view;

    // m_mvp = m_cam_projection * m_mv;
    m_mvp = projection * view;

    // ================================================================

    set_face_uniforms();
    set_pl_uniforms();
    set_clipping_uniforms();
  }

  void Basic_viewer::set_face_uniforms()
  {
    m_face_shader.use();

    m_face_shader.setMatrix4f("mvp_matrix", m_mvp.data());
    m_face_shader.setMatrix4f("mv_matrix", m_mv.data());

    m_face_shader.setVec4f("light_pos", m_lightPosition.data());
    m_face_shader.setVec4f("light_diff", m_diffuse.data());
    m_face_shader.setVec4f("light_spec", m_specular.data());
    m_face_shader.setVec4f("light_amb", m_ambient.data());
    m_face_shader.setFloat("spec_power", m_shininess);

    m_face_shader.setVec4f("clipPlane", m_clip_plane.data());
    m_face_shader.setVec4f("pointPlane", m_point_plane.data());
    m_face_shader.setFloat("rendering_transparency", m_clipping_plane_rendering_transparency);
  }

  void Basic_viewer::set_pl_uniforms()
  {
    m_pl_shader.use();

    m_pl_shader.setVec4f("clipPlane", m_clip_plane.data());
    m_pl_shader.setVec4f("pointPlane", m_point_plane.data());
    m_pl_shader.setMatrix4f("mvp_matrix", m_mvp.data());
    m_pl_shader.setFloat("point_size", m_sizePoints);
  }

  void Basic_viewer::set_clipping_uniforms()
  {
    m_point_plane = m_clipping_matrix * vec4f(0, 0, 0, 1);
    m_clip_plane = m_clipping_matrix * vec4f(0, 0, 1, 0);
    m_plane_shader.use();

    m_plane_shader.setMatrix4f("vp_matrix", m_mvp.data());
    m_plane_shader.setMatrix4f("m_matrix", m_clipping_matrix.data());
  }

  void Basic_viewer::render_scene()
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

    update_uniforms();

    bool half = m_useClippingPlane == CLIPPING_PLANE_SOLID_HALF_ONLY;

    if (m_drawVertices)
    {
      draw_vertices(half ? DRAW_INSIDE_ONLY : DRAW_ALL);
    }
    if (m_drawEdges)
    {
      draw_edges(half ? DRAW_INSIDE_ONLY : DRAW_ALL);
    }
    if (m_drawFaces)
    {
      draw_faces();
    }
    if (m_drawRays)
    {
      draw_rays();
    }
    if (m_drawLines)
    {
      draw_lines();
    }
  }

  vec4f Basic_viewer::color_to_vec4(const CGAL::IO::Color &c) const
  {
    return {(float)c.red() / 255, (float)c.green() / 255, (float)c.blue() / 255, 1.0f};
  }

  void Basic_viewer::draw_faces()
  {
    m_face_shader.use();

    if (m_useClippingPlane == CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF)
    {
      // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
      // Before rendering all transparent objects, disable z-testing first.

      // 1. draw solid first
      draw_faces_(DRAW_INSIDE_ONLY);

      // 2. draw transparent layer second with back face culling to avoid messy triangles
      glDepthMask(false); // disable z-testing
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      glFrontFace(GL_CW);
      draw_faces_(DRAW_OUTSIDE_ONLY);

      // 3. draw solid again without culling and blend to make sure the solid mesh is visible
      glDepthMask(true); // enable z-testing
      glDisable(GL_CULL_FACE);
      glDisable(GL_BLEND);
      draw_faces_(DRAW_INSIDE_ONLY);

      // 4. render clipping plane here
      render_clipping_plane();
      return;
    }

    if (m_useClippingPlane == CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
        m_useClippingPlane == CLIPPING_PLANE_SOLID_HALF_ONLY)
    {
      // 1. draw solid HALF
      draw_faces_(DRAW_INSIDE_ONLY);

      // 2. render clipping plane here
      render_clipping_plane();
      return;
    }

    // 1. draw solid FOR ALL
    draw_faces_(DRAW_ALL);
  }

  void Basic_viewer::draw_faces_(RenderMode mode)
  {
    m_face_shader.use();
    m_face_shader.setFloat("rendering_mode", mode);

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

  void Basic_viewer::draw_rays()
  {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", RenderMode::DRAW_ALL);

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

  void Basic_viewer::draw_vertices(RenderMode render)
  {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", render);

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

  void Basic_viewer::draw_lines()
  {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", RenderMode::DRAW_ALL);

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

  void Basic_viewer::draw_edges(RenderMode mode)
  {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", mode);

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

  void Basic_viewer::generate_clipping_plane()
  {
    size_t size = ((m_scene->bounding_box().xmax() - m_scene->bounding_box().xmin()) +
                   (m_scene->bounding_box().ymax() - m_scene->bounding_box().ymin()) +
                   (m_scene->bounding_box().zmax() - m_scene->bounding_box().zmin()));

    const unsigned int nbSubdivisions = 30;

    auto &array = m_array_for_clipping_plane;
    array.clear();
    for (unsigned int i = 0; i <= nbSubdivisions; ++i)
    {
      const float pos = float(size * (2.0 * i / nbSubdivisions - 1.0));
      array.push_back(pos);
      array.push_back(-float(size));
      array.push_back(0.f);

      array.push_back(pos);
      array.push_back(float(size));
      array.push_back(0.f);

      array.push_back(-float(size));
      array.push_back(pos);
      array.push_back(0.f);

      array.push_back(float(size));
      array.push_back(pos);
      array.push_back(0.f);
    }
  }

  void Basic_viewer::render_clipping_plane()
  {
    if (!m_clipping_plane_rendering || !m_is_opengl_4_3)
      return;
    m_plane_shader.use();
    glBindVertexArray(m_vao[VAO_CLIPPING_PLANE]);
    glLineWidth(0.1f);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_array_for_clipping_plane.size() / 3));
  }

  void Basic_viewer::key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
  {
    Basic_viewer *viewer = static_cast<Basic_viewer *>(glfwGetWindowUserPointer(window));
    viewer->on_key_event(key, scancode, action, mods);
  }

  void Basic_viewer::cursor_callback(GLFWwindow *window, double xpos, double ypo)
  {
    Basic_viewer *viewer = static_cast<Basic_viewer *>(glfwGetWindowUserPointer(window));
    viewer->on_cursor_event(xpos, ypo);
  }

  void Basic_viewer::mouse_btn_callback(GLFWwindow *window, int button, int action, int mods)
  {
    Basic_viewer *viewer = static_cast<Basic_viewer *>(glfwGetWindowUserPointer(window));
    viewer->on_mouse_btn_event(button, action, mods);
  }

  void Basic_viewer::window_size_callback(GLFWwindow *window, int width, int height)
  {
    Basic_viewer *viewer = static_cast<Basic_viewer *>(glfwGetWindowUserPointer(window));

    viewer->m_window_size = {width, height};
    viewer->set_cam_mode(viewer->m_camMode);

    glViewport(0, 0, width, height);
  }

  void Basic_viewer::scroll_callback(GLFWwindow *window, double xoffset, double yoffset)
  {
    Basic_viewer *viewer = static_cast<Basic_viewer *>(glfwGetWindowUserPointer(window));
    viewer->on_scroll_event(xoffset, yoffset);
  }

  void Basic_viewer::start_action(ActionEnum action)
  {
    switch (action)
    {
    case CP_ROTATION:
    case CP_TRANSLATION:
    case CP_TRANS_CAM_DIR:
    case MOUSE_TRANSLATE:
    case MOUSE_ROTATE:
      // glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      break;
    }
  }

  void Basic_viewer::action_event(ActionEnum action)
  {
    if (action == EXIT)
    {
      m_pl_shader.destroy();
      m_face_shader.destroy();
      m_plane_shader.destroy();
      glDeleteBuffers(NB_GL_BUFFERS, m_vbo);
      glDeleteVertexArrays(NB_VAO_BUFFERS, m_vao);
      glfwDestroyWindow(m_window);
      glfwTerminate();
      exit(EXIT_SUCCESS);
    }

    double dt = get_delta_time();
    switch (action)
    {
    case UP:
      // m_camera.translation(0, dt);
      m_camera.move_up(dt);
      break;
    case DOWN:
      // m_camera.translation(0, -dt);
      m_camera.move_down(dt);
      break;
    case LEFT:
      // m_camera.translation(-dt, 0);
      m_camera.move_left(dt);
      break;
    case RIGHT:
      // m_camera.translation(dt, 0);
      m_camera.move_right(dt);
      break;
    case FORWARD:
      m_camera.move(dt);
      break;
    case BACKWARDS:
      m_camera.move(-dt);
      break;
    case MOUSE_ROTATE:
      mouse_rotate();
      break;
    case MOUSE_TRANSLATE:
      mouse_translate();
      break;
    case SWITCH_CAM_MODE:
      m_camera.toggle_mode();
      break;
    case SWITCH_CAM_ROTATION:
      m_camera.toggle_fly();
      // switch_rotation_mode();
      break;
    case FULLSCREEN:
      fullscreen();
      break;
    case SCREENSHOT:
      screenshot("./screenshot.png");
      std::cout << "Screenshot saved in local directory." << std::endl;
      break;
    case INC_MOVE_SPEED_1:
      m_camera.inc_tspeed();
      break;
    case DEC_MOVE_SPEED_1:
      m_camera.dec_tspeed();
      break;
    case INC_ROT_SPEED_1:
      m_camera.inc_rspeed();
      // m_scene_rotation_speed = std::min(0.25f, m_scene_rotation_speed + 0.01f);
      break;
    case DEC_ROT_SPEED_1:
      m_camera.dec_rspeed();
      // m_scene_rotation_speed = std::max(0.25f, m_scene_rotation_speed - 0.01f);
      break;
    case RESET_CAM:
      m_camera.reset_all();
      break;
    case CLIPPING_PLANE_DISPLAY:
      m_clipping_plane_rendering = !m_clipping_plane_rendering;
      break;
    case CLIPPING_PLANE_MODE:
      m_useClippingPlane = static_cast<ClippingMode>((m_useClippingPlane + 1) % CLIPPING_PLANE_END_INDEX);
      break;
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
      m_sizeEdges = std::min(50.f, m_sizeEdges + 0.1f);
      break;
    case DEC_EDGES_SIZE:
      m_sizeEdges = std::max(2.f, m_sizeEdges - 0.1f);
      break;
    case INC_POINTS_SIZE:
      m_sizePoints = std::min(15.f, m_sizePoints + 0.1f);
      break;
    case DEC_POINTS_SIZE:
      m_sizePoints = std::max(5.f, m_sizePoints - 0.1f);
      break;
    case INC_LIGHT_ALL:
      m_ambient.x() += 0.01;
      if (m_ambient.x() > 1)
        m_ambient.x() = 1;
      m_ambient.y() += 0.01;
      if (m_ambient.y() > 1)
        m_ambient.y() = 1;
      m_ambient.z() += 0.01;
      if (m_ambient.z() > 1)
        m_ambient.z() = 1;
      break;
    case DEC_LIGHT_ALL:
      m_ambient.x() -= 0.01;
      if (m_ambient.x() < 0)
        m_ambient.x() = 0;
      m_ambient.y() -= 0.01;
      if (m_ambient.y() < 0)
        m_ambient.y() = 0;
      m_ambient.z() -= 0.01;
      if (m_ambient.z() < 0)
        m_ambient.z() = 0;
      break;
    case INC_LIGHT_R:
      m_ambient.x() += 0.01;
      if (m_ambient.x() > 1)
        m_ambient.x() = 1;
      break;
    case INC_LIGHT_G:
      m_ambient.y() += 0.01;
      if (m_ambient.y() > 1)
        m_ambient.y() = 1;
      break;
    case INC_LIGHT_B:
      m_ambient.z() += 0.01;
      if (m_ambient.z() > 1)
        m_ambient.z() = 1;
      break;
    case DEC_LIGHT_R:
      m_ambient.x() -= 0.01;
      if (m_ambient.x() < 0)
        m_ambient.x() = 0;
      break;
    case DEC_LIGHT_G:
      m_ambient.y() -= 0.01;
      if (m_ambient.y() < 0)
        m_ambient.y() = 0;
      break;
    case DEC_LIGHT_B:
      m_ambient.z() -= 0.01;
      if (m_ambient.z() < 0)
        m_ambient.z() = 0;
      break;
    case CP_ROTATION:
      rotate_clipping_plane();
      break;
    case CP_TRANSLATION:
      translate_clipping_plane();
      break;
    case CP_TRANS_CAM_DIR:
      translate_clipping_plane_cam_dir();
      break;
    case CONSTRAINT_AXIS:
      m_cstr_axis_enum = (m_cstr_axis_enum + 1) % NB_AXIS_ENUM;
      switch_axis(m_cstr_axis_enum);
      break;
    }
  }

  void Basic_viewer::end_action(ActionEnum action)
  {
    switch (action)
    {
    case MOUSE_TRANSLATE:
    case MOUSE_ROTATE:
      glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
      break;
    }
  }

  void Basic_viewer::double_click_event(int btn)
  {
    if (btn == GLFW_MOUSE_BUTTON_2)
    {
      m_camera.reset_position();
    }
  }

  void Basic_viewer::init_keys_actions()
  {
    add_action(GLFW_KEY_ESCAPE, false, EXIT);

    add_action(GLFW_KEY_UP, GLFW_KEY_LEFT_SHIFT, true, FORWARD);
    add_action(GLFW_KEY_DOWN, GLFW_KEY_LEFT_SHIFT, true, BACKWARDS);

    add_action(GLFW_KEY_UP, true, UP);
    add_action(GLFW_KEY_DOWN, true, DOWN);
    add_action(GLFW_KEY_LEFT, true, LEFT);
    add_action(GLFW_KEY_RIGHT, true, RIGHT);

    add_action(GLFW_KEY_O, false, SWITCH_CAM_MODE);
    add_action(GLFW_KEY_V, GLFW_KEY_LEFT_CONTROL, false, SWITCH_CAM_ROTATION);

    add_action(GLFW_KEY_ENTER, GLFW_KEY_LEFT_ALT, false, FULLSCREEN);
    add_action(GLFW_KEY_F1, false, SCREENSHOT);

    add_action(GLFW_KEY_X, false, INC_MOVE_SPEED_1);
    add_action(GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT, false, DEC_MOVE_SPEED_1);

    add_action(GLFW_KEY_R, false, INC_ROT_SPEED_1);
    add_action(GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT, false, DEC_ROT_SPEED_1);

    add_mouse_action(GLFW_MOUSE_BUTTON_1, true, MOUSE_ROTATE);
    add_mouse_action(GLFW_MOUSE_BUTTON_2, true, MOUSE_TRANSLATE);

    add_action(GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, false, RESET_CAM);
    // add_mouse_action(GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, false, CENTER_CAM);

    add_action(GLFW_KEY_C, false, CLIPPING_PLANE_MODE);
    add_action(GLFW_KEY_C, GLFW_KEY_LEFT_ALT, false, CLIPPING_PLANE_DISPLAY);

    add_action(GLFW_KEY_F, false, FACES_DISPLAY);
    add_action(GLFW_KEY_V, false, VERTICES_DISPLAY);
    add_action(GLFW_KEY_E, false, EDGES_DISPLAY);

    add_action(GLFW_KEY_S, false, SHADING_MODE);
    add_action(GLFW_KEY_N, false, INVERSE_NORMAL);
    add_action(GLFW_KEY_M, false, MONO_COLOR);

    add_action(GLFW_KEY_H, GLFW_KEY_LEFT_CONTROL, true, DEC_POINTS_SIZE);
    add_action(GLFW_KEY_J, GLFW_KEY_LEFT_CONTROL, true, INC_POINTS_SIZE);
    add_action(GLFW_KEY_H, true, DEC_EDGES_SIZE);
    add_action(GLFW_KEY_J, true, INC_EDGES_SIZE);

    add_action(GLFW_KEY_PAGE_UP, true, INC_LIGHT_ALL);
    add_action(GLFW_KEY_PAGE_DOWN, true, DEC_LIGHT_ALL);

    add_action(GLFW_KEY_PAGE_UP, GLFW_KEY_LEFT_SHIFT, true, INC_LIGHT_R);
    add_action(GLFW_KEY_PAGE_DOWN, GLFW_KEY_LEFT_SHIFT, true, DEC_LIGHT_R);
    add_action(GLFW_KEY_PAGE_UP, GLFW_KEY_LEFT_ALT, true, INC_LIGHT_G);
    add_action(GLFW_KEY_PAGE_DOWN, GLFW_KEY_LEFT_ALT, true, DEC_LIGHT_G);
    add_action(GLFW_KEY_PAGE_UP, GLFW_KEY_LEFT_CONTROL, true, INC_LIGHT_B);
    add_action(GLFW_KEY_PAGE_DOWN, GLFW_KEY_LEFT_CONTROL, true, DEC_LIGHT_B);

    add_mouse_action(GLFW_MOUSE_BUTTON_1, GLFW_KEY_LEFT_CONTROL, true, CP_ROTATION);
    add_action(GLFW_KEY_A, GLFW_KEY_LEFT_CONTROL, false, CONSTRAINT_AXIS);

    add_mouse_action(GLFW_MOUSE_BUTTON_2, GLFW_KEY_LEFT_CONTROL, true, CP_TRANSLATION);
    add_mouse_action(GLFW_MOUSE_BUTTON_MIDDLE, GLFW_KEY_LEFT_CONTROL, true, CP_TRANS_CAM_DIR);

    /*===================== BIND DESCRIPTIONS ============================*/

    set_action_description({{FORWARD, "Move forward"},
                            {BACKWARDS, "Move backwards"},
                            {UP, "Move right"},
                            {RIGHT, "Move right"},
                            {LEFT, "Move left"},
                            {DOWN, "Move down"},

                            {SWITCH_CAM_MODE, "Switch to Perspective/Orthographic view"},
                            {SWITCH_CAM_ROTATION, "Switch to default/first person mode"},

                            {FULLSCREEN, "Switch to windowed/fullscreen mode"},
                            {SCREENSHOT, "Take a screenshot of the current view"},

                            {MOUSE_ROTATE, "Rotate the view"},
                            {MOUSE_TRANSLATE, "Move the view"},
                            {RESET_CAM, "Reset camera focus"},

                            {CLIPPING_PLANE_MODE, "Switch clipping plane display mode"},
                            {CLIPPING_PLANE_DISPLAY, "Toggle clipping plane rendering on/off"},
                            {CP_ROTATION, "Rotate the clipping plane when enabled"},
                            {CONSTRAINT_AXIS, "Toggle constraint axis for clipping plane rotation"},
                            {CP_TRANSLATION, "Translate the clipping plane when enabled"},
                            {CP_TRANS_CAM_DIR, "Translate the clipping plane along camera direction axis when enabled"},

                            {INC_LIGHT_ALL, "Increase light (all colors, use shift/alt/ctrl for one rgb component)"},
                            {DEC_LIGHT_ALL, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)"},

                            {VERTICES_DISPLAY, "Toggles vertices display"},
                            {EDGES_DISPLAY, "Toggles edges display"},
                            {FACES_DISPLAY, "Toggles faces display"},

                            {INC_POINTS_SIZE, "Increase size of vertices"},
                            {DEC_POINTS_SIZE, "Decrease size of vertices"},
                            {INC_EDGES_SIZE, "Increase size of edges"},
                            {DEC_EDGES_SIZE, "Decrease size of edges"},

                            {MONO_COLOR, "Toggles mono color"},
                            {INVERSE_NORMAL, "Inverse direction of normals"},
                            {SHADING_MODE, "Switch between flat/Gouraud shading display"},
                            {EXIT, "Exits program"}});
  }

  void Basic_viewer::switch_axis(int axis)
  {
    if (axis == X_AXIS)
    {
      std::cout << "Constrained on X" << std::endl;
      m_cstr_axis << 1., 0., 0.;
      return;
    }
    if (axis == Y_AXIS)
    {
      std::cout << "Constrained on Y" << std::endl;
      m_cstr_axis << 0., 1., 0.;
      return;
    }
    if (axis == Z_AXIS)
    {
      std::cout << "Constrained on Z" << std::endl;
      m_cstr_axis << 0., 0., 1.;
      return;
    }
    std::cout << "Constraint Axis Disabled" << std::endl;
  }

  // Normalize Device Coordinates
  vec2f Basic_viewer::to_ndc(double x, double y)
  {
    vec2f result;
    result << x / m_window_size.x() * 2 - 1,
        y / m_window_size.y() * 2 - 1;

    return result;
  }

  // mouse position mapped to the hemisphere
  vec3f Basic_viewer::mapping_cursor_toHemisphere(double x, double y)
  {
    vec3f pt{x, y, 0.};
    float xy_squared = pt.x() * pt.x() + pt.y() * pt.y();
    if (xy_squared > .5)
    { // inside the sphere
      pt.z() = .5 / sqrt(xy_squared);
      pt.normalize();
    }
    else
    {
      // x²+y²+z²=r² => z²=r²-x²-y²
      pt.z() = sqrt(1. - xy_squared);
    }

    if (m_cstr_axis_enum == NO_AXIS)
      return pt;

    // projection on the constraint axis
    float dot = pt.dot(m_cstr_axis);
    vec3f proj = pt - (m_cstr_axis * dot);
    float norm = proj.norm();

    if (norm > 0.)
    {
      float s = 1. / norm;
      if (proj.z() < 0.)
        s = -s;
      pt = proj * s;
    }
    else if (m_cstr_axis.z() == 1.)
    {
      pt << 1., 0., 0.;
    }
    else
    {
      pt << -m_cstr_axis.y(), m_cstr_axis.x(), 0.;
      pt.normalize();
    }

    return pt;
  }

  mat4f Basic_viewer::get_rotation(vec3f const &start, vec3f const &end)
  {
    vec3f rotation_axis = start.cross(end).normalized();
    const float dot = start.dot(end);
    const float angle = acos(std::min(1.f, dot));

    const float d = m_clipping_plane_rot_speed;
    // std::cout << "theta angle : " << angle << std::endl;
    Eigen::Affine3f transform{Eigen::AngleAxisf(angle * d, rotation_axis).toRotationMatrix()};
    return transform.matrix();
  }

  /*********************CLIP STUFF**********************/

  void Basic_viewer::rotate_clipping_plane()
  {
    vec2f cursor_old = get_cursor_old();
    vec2f cursor_current = get_cursor();

    if (cursor_current.x() == cursor_old.x() &&
        cursor_current.y() == cursor_old.y())
      return;

    vec2f old_pos = to_ndc(cursor_old.x(), cursor_old.y());
    vec3f start = mapping_cursor_toHemisphere(old_pos.x(), old_pos.y());

    vec2f crr_pos = to_ndc(cursor_current.x(), cursor_current.y());
    vec3f end = mapping_cursor_toHemisphere(crr_pos.x(), crr_pos.y());

    mat4f rotation = get_rotation(start, end);
    m_clipping_matrix = rotation * m_clipping_matrix;
  }

  void Basic_viewer::translate_clipping_plane()
  {
    vec2f mouse_current = get_cursor();
    const float d = m_clipping_plane_move_speed;

    vec2f delta = get_cursor_delta();
    vec3f dir;
    dir << delta.x(), delta.y(), 0.0f;

    vec3f up{0, 1, 0};
    vec3f right = (-up.cross(m_camForward)).normalized();
    up = m_camForward.cross(right).normalized();

    vec3f result =
        dir.x() * right * d +
        dir.y() * up * d;

    Eigen::Affine3f transform{Eigen::Translation3f(result)};
    mat4f translation = transform.matrix();
    m_clipping_matrix = translation * m_clipping_matrix;
  }

  void Basic_viewer::translate_clipping_plane_cam_dir()
  {

    vec2f cursor_delta = get_cursor_delta();

    float s = cursor_delta.x();
    if (abs(cursor_delta.y()) > abs(cursor_delta.x()))
      s = -cursor_delta.y();

    s *= m_clipping_plane_move_speed;
    Eigen::Affine3f transform{Eigen::Translation3f(s * m_camera.forward())};
    mat4f translation = transform.matrix();
    m_clipping_matrix = translation * m_clipping_matrix;
  }

  /*********************CAM STUFF**********************/

  void Basic_viewer::translate(vec3f dir)
  {
    const float delta = 1.0f / 60;
    vec3f right = vec3f(0, 1, 0).cross(m_camForward).normalized();
    vec3f up = m_camForward.cross(right);

    vec3f result =
        dir.x() * right * delta +
        dir.y() * up * delta +
        dir.z() * m_camForward * delta;

    m_camPosition += result;
  }

  void Basic_viewer::mouse_rotate()
  {
    vec2f delta = get_cursor_delta();
    
    double dt = get_delta_time();

    m_camera.rotation(
      dt * delta.x(), 
      dt * delta.y()
    );
    // std::cout << m_camera.position().x() << " "
    //   << m_camera.position().y() << " "
    //   << m_camera.position().z() << std::endl;
    // std::cout
    //   << m_camera.forward().x() << " "
    //   << m_camera.forward().y() << " "
    //   << m_camera.forward().z() <<
    // std::endl;
  }

  // void Basic_viewer::mouse_rotate(){
  //   vec2f cursor_delta = radians(get_cursor_delta());

  //   if (m_cam_rotation_mode == FREE){
  //     m_cam_view += cursor_delta * m_cam_rotation_speed;

  //     m_camForward = sphericalToCartesian(m_cam_view);

  //     return;
  //   }

  //   m_scene_view += cursor_delta * m_scene_rotation_speed;
  //   m_scene_rotation = eulerAngleXY(-m_scene_view.y(), m_scene_view.x());
  // }

  void Basic_viewer::set_cam_mode(CAM_MODE mode)
  {
    m_camMode = mode;

    float ratio = (float)m_window_size.x() / m_window_size.y();

    if (m_camMode == PERSPECTIVE)
    {
      m_cam_projection = perspective(radians(45.f), ratio, 0.1f, 1000.0f);
      return;
    }

    m_cam_projection = ortho(0.0f, m_cam_orth_zoom * ratio, 0.0f, m_cam_orth_zoom, 0.1f, 100.0f);
  }

  void Basic_viewer::switch_rotation_mode()
  {
    m_cam_rotation_mode = m_cam_rotation_mode == FREE ? OBJECT : FREE;

    if (m_cam_rotation_mode == FREE)
    {
      m_cam_view = cartesianToSpherical(m_camForward);
    }
  }

  void Basic_viewer::fullscreen()
  {
    m_is_fullscreen = !m_is_fullscreen;

    if (m_is_fullscreen)
    {
      int count;
      m_old_window_size = m_window_size;

      GLFWmonitor *monitor = glfwGetMonitors(&count)[0];
      const GLFWvidmode *mode = glfwGetVideoMode(monitor);

      glfwGetWindowPos(m_window, &m_old_window_pos.x(), &m_old_window_pos.y());
      glfwSetWindowMonitor(m_window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
      glViewport(0, 0, mode->width, mode->height);

      std::cout << m_window_size.x() << " " << m_window_size.y();
      return;
    }

    m_window_size = m_old_window_size;
    glfwSetWindowMonitor(m_window, nullptr, m_old_window_pos.x(), m_old_window_pos.y(), m_window_size.x(), m_window_size.y(), 60);
    glViewport(0, 0, m_window_size.x(), m_window_size.y());
  }

  void Basic_viewer::mouse_translate()
  {
    vec2f delta = get_cursor_delta();

    double dt = get_delta_time();

    m_camera.translation(
      dt * -delta.x(),
      dt * delta.y());
  }

  // void Basic_viewer::mouse_translate(){
  //   vec2f delta2 = get_cursor_delta();
  //   vec3f cursor_delta;
  //   cursor_delta << delta2.x(), delta2.y(), 0;

  //   if (cursor_delta.x() == 0 && cursor_delta.y() == 0)
  //     return;

  //   translate(cursor_delta.normalized() * m_cam_speed);
  // }

  void Basic_viewer::print_help()
  {
    std::map<Input::ActionEnum, std::vector<KeyData>> action_keys = get_action_keys();

    std::cout << std::endl
              << "Help for Basic Viewer OpenGl :" << std::endl;

    for (auto pair : action_keys)
    {
      std::vector<KeyData> shortcuts = pair.second;
      ActionEnum action = pair.first;

      std::string line;

      std::string action_str = get_action_description(action);

      // Skip this entry if it has no useful description
      if (action_str.empty())
        continue;

      line += "   " + action_str;

      if (shortcuts.size() > 1)
      {
        line += " (Alternatives : ";

        line += get_key_string(shortcuts[1]) + " ";

        for (int s = 2; s < shortcuts.size(); s++)
        {
          line += ", " + get_key_string(shortcuts[s]);
        }

        line += ").";
      }

      std::cout
          << std::setw(12)
          << (shortcuts.size() > 0 ? get_key_string(shortcuts[0]) : "(unbound)")
          << line
          << std::setw(0)
          << std::endl;
    }
  }

  void Basic_viewer::scroll_event()
  {
    double yoffset = get_scroll_yoffset();
    
    int state = glfwGetKey(m_window, GLFW_KEY_LEFT_CONTROL);

    if (state == GLFW_PRESS) {
      m_camera.modify_fov(yoffset);
    } else {
      m_camera.move(8.f * yoffset * get_delta_time());
    }
  }

  void Basic_viewer::screenshot(const std::string &filepath)
  {
    // https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/ (thanks)
    // https://github.com/nothings/stb/
    // The stb lib used here is from glfw/deps

    const GLsizei nrChannels = 3;
    GLsizei stride = nrChannels * m_window_size.x();
    stride += (stride % 4) ? (4 - stride % 4) : 0; // stride must be a multiple of 4
    GLsizei bufferSize = stride * 3 * m_window_size.y();

    std::vector<char> buffer(bufferSize);

    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, m_window_size.x(), m_window_size.y(), GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath.data(), m_window_size.x(), m_window_size.y(), nrChannels, buffer.data(), stride);
  }
} // end namespace GLFW 
} // end namespace CGAL
