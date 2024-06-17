#pragma once

#include "Basic_viewer.h"

namespace CGAL::GLFW {
  void Basic_Viewer::error_callback(int error, const char *description)
  {
    //fprintf(stderr, "GLFW returned an error:\n\t%s (%i)\n", description, error);
  }

  GLFWwindow* Basic_Viewer::create_window(int width, int height, const char *title, bool hidden)
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


  Basic_Viewer::Basic_Viewer(
      const Graphics_scene* graphics_scene,
      const char *title,
      bool draw_vertices,
      bool draw_edges,
      bool draw_faces,
      bool use_mono_color,
      bool inverse_normal,
      bool draw_rays,
      bool draw_lines) : 
    m_scene(graphics_scene),
    m_title(title),
    m_draw_vertices(draw_vertices),
    m_draw_edges(draw_edges),
    m_draw_rays(draw_rays),
    m_draw_lines(draw_lines),
    m_draw_faces(draw_faces),
    m_use_mono_color(use_mono_color),
    m_inverse_normal(inverse_normal)
    {
      init_keys_actions();
    }

    void Basic_Viewer::show()
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
      set_cam_mode(m_cam_mode);

      GLint major, minor;
      glGetIntegerv(GL_MAJOR_VERSION, &major);
      glGetIntegerv(GL_MINOR_VERSION, &minor);

      if (major > 4 || major == 4 && minor >= 3){
        m_is_opengl_4_3 = true;
      }

      compile_shaders();

      vec3f pmin(
        m_scene->bounding_box().xmin(), 
        m_scene->bounding_box().ymin(),
        m_scene->bounding_box().zmin()
      );

      vec3f pmax(
        m_scene->bounding_box().xmax(), 
        m_scene->bounding_box().ymax(),
        m_scene->bounding_box().zmax()
      );
      
      m_camera.lookat(pmin, pmax);

      while (!glfwWindowShouldClose(m_window))
      {
        render_scene();
        glfwSwapBuffers(m_window);
        handle_events();
      }

      glfwTerminate();
    }

    void Basic_Viewer::make_screenshot(const std::string& pngpath) {
      m_are_buffers_initialized = false;
      m_window = create_window(m_window_size.x(), m_window_size.y(), m_title, true);
      init_buffers();

      set_cam_mode(m_cam_mode); 
      
      GLint major, minor;
      glGetIntegerv(GL_MAJOR_VERSION, &major);
      glGetIntegerv(GL_MINOR_VERSION, &minor);

      if (major > 4 || major == 4 && minor >= 3){
        m_is_opengl_4_3 = true;
      }

      compile_shaders();
      render_scene();
      glfwSwapBuffers(m_window);
      screenshot(pngpath);
      glfwTerminate();
    }

  void Basic_Viewer::compile_shaders() { 
    const char* face_vert = m_is_opengl_4_3 ? vertex_source_color : vertex_source_color_comp;
    const char* face_frag = m_is_opengl_4_3 ? fragment_source_color : fragment_source_color_comp;
    const char* pl_vert = m_is_opengl_4_3 ? vertex_source_p_l : vertex_source_p_l_comp;
    const char* pl_frag = m_is_opengl_4_3 ? fragment_source_p_l : fragment_source_p_l_comp;
    const char* plane_vert = vertex_source_clipping_plane;
    const char* plane_frag = fragment_source_clipping_plane;


    m_face_shader = Shader::loadShader(face_vert, face_frag, "FACE");
    m_pl_shader = Shader::loadShader(pl_vert, pl_frag, "PL");
    m_plane_shader = Shader::loadShader(plane_vert, plane_frag, "PLANE");
  }

  void Basic_Viewer::load_buffer(int i, int location, const std::vector<float>& vector, int dataCount){
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[i]);

    glBufferData(GL_ARRAY_BUFFER, vector.size() * sizeof(float), vector.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(location, dataCount, GL_FLOAT, GL_FALSE, dataCount * sizeof(float), nullptr);

    glEnableVertexAttribArray(location);
  }


  void Basic_Viewer::load_buffer(int i, int location, int gsEnum, int dataCount){ 
    const std::vector<float>& vector = m_scene->get_array_of_index(gsEnum);
    load_buffer(i, location, vector, dataCount);
  }

  void Basic_Viewer::init_buffers(){
    if (m_are_buffers_initialized){
      glGenBuffers(NB_GL_BUFFERS, m_vbo);
      glGenVertexArrays(NB_VAO_BUFFERS, m_vao); 
      m_are_buffers_initialized = true;
    }
  }

  void Basic_Viewer::load_scene()
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
    if (m_flat_shading) {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_MONO_FACES, 3);
    } else {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_MONO_FACES, 3);
    }

    // 5.2) Colored faces
    glBindVertexArray(m_vao[VAO_COLORED_FACES]); 
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_FACES, 3);
    if (m_flat_shading) {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_COLORED_FACES, 3);
    } else {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_COLORED_FACES, 3);
    }
    load_buffer(bufn++, 2, Graphics_scene::COLOR_FACES, 3);

    // 6) clipping plane shader
    if (m_is_opengl_4_3) {
      generate_clipping_plane();
      glBindVertexArray(m_vao[VAO_CLIPPING_PLANE]);
      load_buffer(bufn++, 0, m_array_for_clipping_plane, 3);
    }

    m_are_buffers_initialized = true;
  }

  CGAL::Plane_3<Basic_Viewer::Local_kernel> Basic_Viewer::clipping_plane() const
  {
    const mat4f cpm = m_clipping_matrix;
    CGAL::Aff_transformation_3<Basic_Viewer::Local_kernel> aff(
      cpm(0,0), cpm(0,1), cpm(0,2), cpm(0,3),
      cpm(1,0), cpm(1,1), cpm(1,2), cpm(1,3),
      cpm(2,0), cpm(2,1), cpm(2,2), cpm(2,3)
    );

    CGAL::Plane_3<Local_kernel> p3(0, 0, 1, 0);
    return p3.transform(aff);
  }

  void Basic_Viewer::update_uniforms(){
    m_model_view = lookAt(m_cam_position, m_cam_position + m_cam_forward, vec3f(0,1,0)) * m_scene_rotation;

    mat4f view = m_camera.view();
    mat4f projection = m_camera.projection(m_window_size.x(), m_window_size.y(), 45.f);
    
    m_mvp = m_cam_projection * m_model_view;
    // m_mvp =projection * view; 

    // ================================================================

    set_face_uniforms();
    set_pl_uniforms();
    set_clipping_uniforms();
  }

  void Basic_Viewer::set_face_uniforms() {
    m_face_shader.use();

    m_face_shader.setMatrix4f("mvp_matrix", m_mvp.data());
    m_face_shader.setMatrix4f("mv_matrix", m_model_view.data());
    
    m_face_shader.setVec4f("light_pos", m_light_position.data());
    m_face_shader.setVec4f("light_diff", m_diffuse.data());
    m_face_shader.setVec4f("light_spec", m_specular.data());
    m_face_shader.setVec4f("light_amb", m_ambient.data());
    m_face_shader.setFloat("spec_power", m_shininess);    

    m_face_shader.setVec4f("clipPlane", m_clip_plane.data());
    m_face_shader.setVec4f("pointPlane", m_point_plane.data());
    m_face_shader.setFloat("rendering_transparency", m_clipping_plane_rendering_transparency);
  }

  void Basic_Viewer::set_pl_uniforms() {
    m_pl_shader.use();
    
    m_pl_shader.setVec4f("clipPlane", m_clip_plane.data());
    m_pl_shader.setVec4f("pointPlane", m_point_plane.data());
    m_pl_shader.setMatrix4f("mvp_matrix", m_mvp.data());
    m_pl_shader.setFloat("point_size", m_size_points);
  }

  void Basic_Viewer::set_clipping_uniforms() {
    m_point_plane = m_clipping_matrix * vec4f(0, 0, 0, 1);
    m_clip_plane = m_clipping_matrix * vec4f(0, 0, 1, 0);
    m_plane_shader.use();

    m_plane_shader.setMatrix4f("vp_matrix", m_mvp.data());
    m_plane_shader.setMatrix4f("m_matrix", m_clipping_matrix.data());
  }
  
  void Basic_Viewer::render_scene()
  {
    if(!m_are_buffers_initialized) { load_scene(); }
    
    glClearColor(1.0f,1.0f,1.0f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);
    
    update_uniforms();

    bool half = m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY;
    
    if (m_draw_vertices)  { 
      draw_vertices(half ? DRAW_INSIDE_ONLY : DRAW_ALL); 
    }
    if (m_draw_edges) {
      draw_edges(half ? DRAW_INSIDE_ONLY : DRAW_ALL); 
    }
    if (m_draw_faces)     { draw_faces(); }
    if (m_draw_rays)      { draw_rays(); } 
    if (m_draw_lines)     { draw_lines(); }
  }

  vec4f Basic_Viewer::color_to_vec4(const CGAL::IO::Color& c) const
  {
    return { (float)c.red()/255, (float)c.green()/255, (float)c.blue()/255, 1.0f };
  }

  void Basic_Viewer::draw_faces() 
  {
    m_face_shader.use();

    if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF) {
      // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
      // Before rendering all transparent objects, disable z-testing first.

      // 1. draw solid first
      draw_faces_(DRAW_INSIDE_ONLY);

      // 2. draw transparent layer second with back face culling to avoid messy triangles
      glDepthMask(false); //disable z-testing
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      glFrontFace(GL_CW);
      draw_faces_(DRAW_OUTSIDE_ONLY);

      // 3. draw solid again without culling and blend to make sure the solid mesh is visible
      glDepthMask(true); //enable z-testing
      glDisable(GL_CULL_FACE);
      glDisable(GL_BLEND);
      draw_faces_(DRAW_INSIDE_ONLY);

      // 4. render clipping plane here
      render_clipping_plane();
      return;
    }
      
    if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
        m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) 
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

  void Basic_Viewer::draw_faces_(RenderMode mode){
    m_face_shader.use();
    m_face_shader.setFloat("rendering_mode", mode);

    vec4f color = color_to_vec4(m_faces_mono_color);    

    glBindVertexArray(m_vao[VAO_MONO_FACES]);
    glVertexAttrib4fv(2, color.data());
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_FACES));
  
    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    if (m_use_mono_color) {
      glDisableVertexAttribArray(2);
    } else {
      glEnableVertexAttribArray(2);
    }

    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_FACES));
  }

  void Basic_Viewer::draw_rays() {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", RenderMode::DRAW_ALL);
    
    vec4f color = color_to_vec4(m_rays_mono_color);    

    glBindVertexArray(m_vao[VAO_MONO_RAYS]);
    glVertexAttrib4fv(1, color.data());

    glLineWidth(m_size_rays);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_RAYS));
  
    glBindVertexArray(m_vao[VAO_COLORED_RAYS]);
    if (m_use_mono_color) {
      glDisableVertexAttribArray(1);
    } else {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_RAYS));
  }

  void Basic_Viewer::draw_vertices(RenderMode render) {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", render);
    
    vec4f color = color_to_vec4(m_vertices_mono_color);    

    glBindVertexArray(m_vao[VAO_MONO_POINTS]);
    glVertexAttrib4fv(1, color.data());
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_POINTS));
  
    glBindVertexArray(m_vao[VAO_COLORED_POINTS]);
    if (m_use_mono_color) {
      glDisableVertexAttribArray(1);
    } else {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_POINTS));

  }

  void Basic_Viewer::draw_lines() {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", RenderMode::DRAW_ALL);
    
    vec4f color = color_to_vec4(m_lines_mono_color);    
    
    glBindVertexArray(m_vao[VAO_MONO_LINES]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_size_lines);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_LINES));
  
  
    glBindVertexArray(m_vao[VAO_COLORED_LINES]);
    if (m_use_mono_color) {
      glDisableVertexAttribArray(1);
    } else {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_LINES));
  }

  void Basic_Viewer::draw_edges(RenderMode mode) {
    m_pl_shader.use();
    m_pl_shader.setFloat("rendering_mode", mode);
          
    vec4f color = color_to_vec4(m_edges_mono_color);    

    glBindVertexArray(m_vao[VAO_MONO_SEGMENTS]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_size_edges);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_SEGMENTS));
  
    glBindVertexArray(m_vao[VAO_COLORED_SEGMENTS]);
    if (m_use_mono_color) {
      glDisableVertexAttribArray(1);
    } else {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_SEGMENTS));
    
  }

  void Basic_Viewer::generate_clipping_plane() {
      size_t size=((m_scene->bounding_box().xmax()-m_scene->bounding_box().xmin()) +
                (m_scene->bounding_box().ymax()-m_scene->bounding_box().ymin()) +
                (m_scene->bounding_box().zmax()-m_scene->bounding_box().zmin()));
      
      const unsigned int nbSubdivisions=30;

      auto& array = m_array_for_clipping_plane;
      array.clear();
      for (unsigned int i=0; i<=nbSubdivisions; ++i)
      {
        const float pos = float(size*(2.0*i/nbSubdivisions-1.0));
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

  void Basic_Viewer::render_clipping_plane() {
    if (!m_clipping_plane_rendering || !m_is_opengl_4_3) return;
    m_plane_shader.use();
    glBindVertexArray(m_vao[VAO_CLIPPING_PLANE]);
    glLineWidth(0.1f);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_array_for_clipping_plane.size()/3));
  }

  void Basic_Viewer::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    Basic_Viewer* viewer = static_cast<Basic_Viewer*>(glfwGetWindowUserPointer(window)); 
    viewer->on_key_event(key, scancode, action, mods);
  }
  
  void Basic_Viewer::cursor_callback(GLFWwindow* window, double xpos, double ypo)
  {
    Basic_Viewer* viewer = static_cast<Basic_Viewer*>(glfwGetWindowUserPointer(window)); 
    viewer->on_cursor_event(xpos, ypo);
  }
  
  void Basic_Viewer::mouse_btn_callback(GLFWwindow* window, int button, int action, int mods)
  {
    Basic_Viewer* viewer = static_cast<Basic_Viewer*>(glfwGetWindowUserPointer(window)); 
    viewer->on_mouse_btn_event(button, action, mods);
  }

  void Basic_Viewer::window_size_callback(GLFWwindow* window, int width, int height) {
    Basic_Viewer* viewer = static_cast<Basic_Viewer*>(glfwGetWindowUserPointer(window)); 

    viewer->m_window_size = {width, height};
    viewer->set_cam_mode(viewer->m_cam_mode);

    glViewport(0, 0, width, height);
  }

  void Basic_Viewer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    Basic_Viewer* viewer = static_cast<Basic_Viewer*>(glfwGetWindowUserPointer(window)); 
    viewer->on_scroll_event(xoffset, yoffset);
  }

  void Basic_Viewer::start_action(ActionEnum action){
    switch (action) {
      case CP_ROTATION:
      case CP_TRANSLATION:
      case CP_TRANS_CAM_DIR:
      case MOUSE_TRANSLATE:
      case MOUSE_ROTATE:
        // glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        break;
    }
  }

  void Basic_Viewer::action_event(ActionEnum action){
    if (action == EXIT) {
      m_pl_shader.destroy();
      m_face_shader.destroy(); 
      m_plane_shader.destroy();
      glDeleteBuffers(NB_GL_BUFFERS, m_vbo);
      glDeleteVertexArrays(NB_VAO_BUFFERS, m_vao);
      glfwDestroyWindow(m_window);
      glfwTerminate();
      exit(EXIT_SUCCESS);
    }

    switch (action){
      case UP:
        translate(vec3f(0, m_cam_speed, 0));
        break;
      case DOWN:
        translate(vec3f(0, -m_cam_speed, 0));
        break;
      case LEFT:
        translate(vec3f(m_cam_speed, 0, 0));
        break;
      case RIGHT:
        translate(vec3f(-m_cam_speed, 0, 0));
        break;
      case FORWARD:
        translate(vec3f(0, 0, m_cam_speed));
        break;
      case BACKWARDS:
        translate(vec3f(0, 0, -m_cam_speed));
        break;
      case MOUSE_ROTATE: 
        mouse_rotate();
        break;
      case MOUSE_TRANSLATE:
        mouse_translate();
        break;
      case SWITCH_CAM_MODE:
        set_cam_mode(m_cam_mode == PERSPECTIVE ? ORTHOGRAPHIC : PERSPECTIVE);
        break;
      case SWITCH_CAM_ROTATION:
        switch_rotation_mode();
        break;
      case FULLSCREEN:
        fullscreen();
        break;
      case SCREENSHOT:
        screenshot("./screenshot.png");
        std::cout << "Screenshot saved in local directory." << std::endl; 
        break;
      case INC_ZOOM:
        zoom(1.0f);
        break;
      case DEC_ZOOM:
        zoom(-1.0f);
        break;
      case INC_MOVE_SPEED_D1:
        m_cam_speed += 0.1f;
        break;
      case DEC_MOVE_SPEED_D1:
        m_cam_speed -= 0.1f;
        break;
      case INC_MOVE_SPEED_1:
        m_cam_speed++;
        break;
      case DEC_MOVE_SPEED_1:
        m_cam_speed--;
        break;
      case INC_ROT_SPEED_D1:
        m_scene_rotation_speed += 0.1f;
        break;
      case DEC_ROT_SPEED_D1:
        m_scene_rotation_speed -= 0.1f;
        break;
      case INC_ROT_SPEED_1:
        m_scene_rotation_speed++;
        break;
      case DEC_ROT_SPEED_1:
        m_scene_rotation_speed--;
        break;
      case CLIPPING_PLANE_DISPLAY:
        m_clipping_plane_rendering = !m_clipping_plane_rendering;
        break;
      case CLIPPING_PLANE_MODE:
        m_use_clipping_plane = static_cast<ClippingMode>((m_use_clipping_plane + 1) % CLIPPING_PLANE_END_INDEX);
        break;
      case VERTICES_DISPLAY:
        m_draw_vertices = !m_draw_vertices;
        break;
      case FACES_DISPLAY:
        m_draw_faces = !m_draw_faces;
        break;
      case EDGES_DISPLAY:
        m_draw_edges = !m_draw_edges;
        break;
      case SHADING_MODE:
        m_flat_shading = !m_flat_shading;
        m_are_buffers_initialized = false;
        break;
      case INVERSE_NORMAL:
        m_inverse_normal = !m_inverse_normal;
        m_scene->reverse_all_normals();
        m_are_buffers_initialized = false;
        break;
      case MONO_COLOR:
        m_use_mono_color = !m_use_mono_color;
        break;
      case INC_EDGES_SIZE: 
        if (m_size_edges < 100)
          m_size_edges += 0.1f;
        break;
      case DEC_EDGES_SIZE:
        if (m_size_edges>1)
          m_size_edges -= 0.1f; 
        break;
      case INC_POINTS_SIZE: 
        if (m_size_points < 100)
          m_size_points += 0.1f;
        break;
      case DEC_POINTS_SIZE:
        if (m_size_points>1)
          m_size_points -= 0.1f; 
        break;
      case INC_LIGHT_ALL:
        m_ambient.x() += 0.01;
        if (m_ambient.x() > 1) m_ambient.x()=1; 
        m_ambient.y() += 0.01;
        if (m_ambient.y() > 1) m_ambient.y()=1; 
        m_ambient.z() += 0.01;
        if (m_ambient.z() > 1) m_ambient.z()=1; 
        break;
      case DEC_LIGHT_ALL:
        m_ambient.x() -= 0.01;
        if (m_ambient.x() < 0) m_ambient.x()=0; 
        m_ambient.y()-= 0.01;
        if (m_ambient.y() < 0) m_ambient.y()=0; 
        m_ambient.z()-= 0.01;
        if (m_ambient.z() < 0) m_ambient.z()=0; 
        break;
      case INC_LIGHT_R:
        m_ambient.x()+= 0.01;
        if (m_ambient.x() > 1) m_ambient.x()=1; 
        break;
      case INC_LIGHT_G:
        m_ambient.y()+= 0.01;
        if (m_ambient.y() > 1) m_ambient.y()=1; 
        break;
      case INC_LIGHT_B:
        m_ambient.z()+= 0.01;
        if (m_ambient.z() > 1) m_ambient.z()=1; 
        break;
      case DEC_LIGHT_R:
        m_ambient.x()-= 0.01;
        if (m_ambient.x() < 0) m_ambient.x()=0; 
        break;
      case DEC_LIGHT_G:  
        m_ambient.y()-= 0.01;
        if (m_ambient.y() < 0) m_ambient.y()=0; 
        break;
      case DEC_LIGHT_B:
        m_ambient.z()-= 0.01;
        if (m_ambient.z() < 0) m_ambient.z()=0;
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
        m_cstr_axis_enum = (m_cstr_axis_enum+1) % NB_AXIS_ENUM;
        switch_axis(m_cstr_axis_enum);
        break;
    }
  }

  void Basic_Viewer::end_action(ActionEnum action){
    switch (action) {
      case MOUSE_TRANSLATE:
      case MOUSE_ROTATE:
        glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        break;
    }
  }

  void Basic_Viewer::init_keys_actions() {
    add_action(GLFW_KEY_ESCAPE, false, EXIT);

    add_action(GLFW_KEY_UP, GLFW_KEY_LEFT_SHIFT, true, FORWARD);
    add_action(GLFW_KEY_DOWN, GLFW_KEY_LEFT_SHIFT, true, BACKWARDS);

    add_action(GLFW_KEY_UP, true, UP);
    add_action(GLFW_KEY_DOWN, true, DOWN);
    add_action(GLFW_KEY_LEFT, true, LEFT);
    add_action(GLFW_KEY_RIGHT, true, RIGHT);

    add_action(GLFW_KEY_O, false, SWITCH_CAM_MODE);
    add_action(GLFW_KEY_V, GLFW_KEY_LEFT_CONTROL, false, SWITCH_CAM_ROTATION);

    add_action(GLFW_KEY_Z, false, INC_ZOOM);
    add_action(GLFW_KEY_Z, GLFW_KEY_LEFT_SHIFT, false, DEC_ZOOM);

    add_action(GLFW_KEY_ENTER, GLFW_KEY_LEFT_ALT, false, FULLSCREEN);
    add_action(GLFW_KEY_F1, false, SCREENSHOT);

    add_action(GLFW_KEY_X, false, INC_MOVE_SPEED_1);
    add_action(GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL, false, INC_MOVE_SPEED_D1);
    add_action(GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT, false, DEC_MOVE_SPEED_1);
    add_action(GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT, GLFW_KEY_LEFT_CONTROL, false, DEC_MOVE_SPEED_D1);

    add_action(GLFW_KEY_R, false, INC_ROT_SPEED_1);
    add_action(GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, false, INC_ROT_SPEED_D1);
    add_action(GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT, false, DEC_ROT_SPEED_1);
    add_action(GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT, GLFW_KEY_LEFT_CONTROL, false, DEC_ROT_SPEED_D1);

    add_mouse_action(GLFW_MOUSE_BUTTON_1, true, MOUSE_ROTATE);
    add_mouse_action(GLFW_MOUSE_BUTTON_2, true, MOUSE_TRANSLATE);
  
    add_action(GLFW_KEY_C, false, CLIPPING_PLANE_MODE);
    add_action(GLFW_KEY_C, GLFW_KEY_LEFT_ALT, false, CLIPPING_PLANE_DISPLAY); 

    add_action(GLFW_KEY_F, false, FACES_DISPLAY);
    add_action(GLFW_KEY_V, false, VERTICES_DISPLAY);
    add_action(GLFW_KEY_E, false, EDGES_DISPLAY);
    // add_action(GLFW_KEY_T, false, TEXT_DISPLAY);
    
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

    set_action_description({
      {FORWARD, "Move forward"},
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
      {EXIT, "Exits program"}
    });
  }
  
  void Basic_Viewer::switch_axis(int axis) {
    if (axis == X_AXIS) {
      std::cout << "Constrained on X" << std::endl;
      m_cstr_axis << 1., 0., 0.;
      return;
    }
    if (axis == Y_AXIS) {
      std::cout << "Constrained on Y" << std::endl;
      m_cstr_axis << 0., 1., 0.;
      return;
    }
    if (axis == Z_AXIS) {
      std::cout << "Constrained on Z" << std::endl;
      m_cstr_axis << 0., 0., 1.;
      return;
    }
    std::cout << "Constraint Axis Disabled" << std::endl;
  } 

  // Normalize Device Coordinates 
  vec2f Basic_Viewer::to_ndc(double x, double y) {
    vec2f result;
    result << 
      x / m_window_size.x() * 2 - 1,
      y / m_window_size.y() * 2 - 1;
    
    return result;
  }

  // mouse position mapped to the hemisphere 
  vec3f Basic_Viewer::mapping_cursor_toHemisphere(double x, double y) {
    vec3f pt { x, y, 0. };
    float xy_squared = pt.x()*pt.x()+pt.y()*pt.y();
    if (xy_squared > .5) { // inside the sphere
      pt.z() = .5/sqrt(xy_squared);
      pt.normalize();
    } else {
      // x²+y²+z²=r² => z²=r²-x²-y²
      pt.z() = sqrt(1. - xy_squared);
    } 

    if (m_cstr_axis_enum == NO_AXIS) return pt;

    // projection on the constraint axis 
    float dot = pt.dot(m_cstr_axis);
    vec3f proj = pt - (m_cstr_axis * dot);
    float norm = proj.norm();

    if (norm > 0.) {
      float s = 1./norm;
      if (proj.z() < 0.) s = -s;
      pt = proj * s; 
    } else if (m_cstr_axis.z() == 1.) {
      pt << 1., 0., 0.;
    } else {
      pt << -m_cstr_axis.y(), m_cstr_axis.x(), 0.;
      pt.normalize();
    }

    return pt;  
  }

  mat4f Basic_Viewer::get_rotation(vec3f const& start, vec3f const& end) {
    vec3f rotation_axis = start.cross(end).normalized();
    const float dot = start.dot(end);
    const float angle = acos(std::min(1.f, dot));

    const float d = m_clipping_plane_rot_speed;
    // std::cout << "theta angle : " << angle << std::endl;
    Eigen::Affine3f transform{Eigen::AngleAxisf(angle*d, rotation_axis).toRotationMatrix()};
    return transform.matrix();
  }

  /*********************CLIP STUFF**********************/

  void Basic_Viewer::rotate_clipping_plane() {
    vec2f cursor_old = get_cursor_old();
    vec2f cursor_current = get_cursor();

    if (cursor_current.x() == cursor_old.x() && 
        cursor_current.y() == cursor_old.y()) return;
  
    vec2f old_pos = to_ndc(cursor_old.x(), cursor_old.y()); 
    vec3f start = mapping_cursor_toHemisphere(old_pos.x(), old_pos.y());

    vec2f crr_pos = to_ndc(cursor_current.x(), cursor_current.y()); 
    vec3f end = mapping_cursor_toHemisphere(crr_pos.x(), crr_pos.y());

    mat4f rotation = get_rotation(start, end);
    m_clipping_matrix = rotation * m_clipping_matrix;
  }

  void Basic_Viewer::translate_clipping_plane() {
    vec2f mouse_current = get_cursor(); 
    const float d = m_clipping_plane_move_speed;

    vec2f delta = get_cursor_delta();
    vec3f dir;
    dir << delta.x(), delta.y(), 0.0f;

    vec3f up {0, 1, 0};
    vec3f right = (-up.cross(m_cam_forward)).normalized(); 
    up = m_cam_forward.cross(right).normalized();

    vec3f result = 
      dir.x() * right * d + 
      dir.y() * up * d;

    Eigen::Affine3f transform { Eigen::Translation3f(result) };
    mat4f translation = transform.matrix(); 
    m_clipping_matrix = translation * m_clipping_matrix;  
  }

  void Basic_Viewer::translate_clipping_plane_cam_dir() {

    vec2f cursor_delta = get_cursor_delta();

    float s = cursor_delta.x();
    if (abs(cursor_delta.y()) > abs(cursor_delta.x()))
      s = -cursor_delta.y();
    
    s *= m_clipping_plane_move_speed;
    Eigen::Affine3f transform { Eigen::Translation3f(s * m_cam_forward) };
    mat4f translation = transform.matrix(); 
    m_clipping_matrix = translation * m_clipping_matrix;  
  }

  /*********************CAM STUFF**********************/

  void Basic_Viewer::translate(vec3f dir){
    const float delta = 1.0f/60;
    vec3f right = vec3f(0, 1, 0).cross(m_cam_forward).normalized(); 
    vec3f up = m_cam_forward.cross(right);

    vec3f result = 
      dir.x() * right * delta +
      dir.y() * up * delta +
      dir.z() * m_cam_forward * delta;

    m_cam_position += result;
  }

  void Basic_Viewer::mouse_rotate(){
    vec2f cursor_delta = radians(get_cursor_delta());

    if (m_cam_rotation_mode == FREE){
      m_cam_view += cursor_delta * m_cam_rotation_speed;

      m_cam_forward = sphericalToCartesian(m_cam_view);

      return;
    }

    m_scene_view += cursor_delta * m_scene_rotation_speed;
    m_scene_rotation = eulerAngleXY(-m_scene_view.y(), m_scene_view.x());
  }

  void Basic_Viewer::set_cam_mode(CAM_MODE mode) {
    m_cam_mode = mode;

    float ratio = (float)m_window_size.x()/m_window_size.y();

    if (m_cam_mode == PERSPECTIVE){
      m_cam_projection = perspective(radians(45.f), ratio, 0.1f, 1000.0f);
      return;
    }
    
    m_cam_projection = ortho(0.0f, m_cam_orth_zoom * ratio, 0.0f, m_cam_orth_zoom, 0.1f, 100.0f);
  }

  void Basic_Viewer::switch_rotation_mode() {
    m_cam_rotation_mode = m_cam_rotation_mode == FREE ? OBJECT : FREE;

    if (m_cam_rotation_mode == FREE) {
      m_cam_view = cartesianToSpherical(m_cam_forward);
    }
  }

  void Basic_Viewer::fullscreen(){
    m_is_fullscreen = !m_is_fullscreen;

    if (m_is_fullscreen) {
      int count;
      m_old_window_size = m_window_size;

      GLFWmonitor* monitor = glfwGetMonitors(&count)[0];
      const GLFWvidmode* mode = glfwGetVideoMode(monitor);

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

  void Basic_Viewer::mouse_translate(){
    vec2f delta2 = get_cursor_delta();
    vec3f cursor_delta;
    cursor_delta << delta2.x(), delta2.y(), 0;

    if (cursor_delta.x() == 0 && cursor_delta.y() == 0)
      return;
    
    translate(cursor_delta.normalized() * m_cam_speed);
  }

  void Basic_Viewer::print_help(){
    std::map<Input::ActionEnum, std::vector<KeyData>> action_keys = get_action_keys();

    std::cout << std::endl << "Help for Basic Viewer OpenGl :" << std::endl;

    for (auto pair : action_keys){
      std::vector<KeyData> shortcuts = pair.second;
      ActionEnum action = pair.first;

      std::string line;

      std::string action_str = get_action_description(action);


      // Skip this entry if it has no useful description
      if (action_str.empty()) continue;
        
      line += "   " + action_str;
      

      if (shortcuts.size() > 1) {
        line += " (Alternatives : ";
        
        line += get_key_string(shortcuts[1]) + " ";

        for (int s = 2; s < shortcuts.size(); s++) {
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

  void Basic_Viewer::zoom(float z){
    if (m_cam_mode == ORTHOGRAPHIC){
      m_cam_orth_zoom += z;
      set_cam_mode(ORTHOGRAPHIC);
      return;
    }

    // readjust position of the camera
    const float zoom_inc = 0.1f; // todo -> define

    m_cam_position += (zoom_inc * z) * m_cam_forward;
    set_cam_mode(PERSPECTIVE);
  }

  void Basic_Viewer::screenshot(const std::string& filepath) {
    // https://lencerf.y()ithub.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/ (thanks)
    // https://github.com/nothings/stb/
    // The stb lib used here is from glfw/deps 
    
    const GLsizei nrChannels = 3;
    GLsizei stride = nrChannels * m_window_size.x();
    stride += (stride % 4) ? (4 - stride % 4) : 0; // stride must be a multiple of 4
    GLsizei bufferSize = stride * m_window_size.y();

    std::vector<char> buffer(bufferSize);

    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, m_window_size.x(), m_window_size.y(), GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath.data(), m_window_size.x(), m_window_size.y(), nrChannels, buffer.data(), stride);
  }

  // Blocking call
  inline void draw_graphics_scene(const Graphics_scene &graphics_scene, const char *title)
  {
    Basic_Viewer(&graphics_scene, title).show();
  }

  inline void draw_graphics_scene(const Graphics_scene *graphics_scene, const char *title)
  {
    Basic_Viewer(graphics_scene, title).show();
  }
} 
