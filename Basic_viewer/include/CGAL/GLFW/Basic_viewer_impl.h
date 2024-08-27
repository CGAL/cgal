#ifndef CGAL_BASIC_VIEWER_GLFW_IMPL_H
#define CGAL_BASIC_VIEWER_GLFW_IMPL_H

#include "Basic_viewer.h"

namespace CGAL 
{
namespace GLFW 
{

  inline 
  void Basic_viewer::error_callback(int error, const char *description)
  {
    std::cerr << "GLFW returned an error:\n\t" << description << "(" << error << ")\n";
  }

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
    if (!hidden)
    {
      glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    }

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

  Basic_viewer::Basic_viewer(
    const Graphics_scene& graphicScene,
    const char* title,
    bool drawVertices,
    bool drawEdges,
    bool drawFaces,
    bool drawRays,
    bool drawLines, 
    bool useDefaultColor,
    bool inverseNormal,
    bool flatShading
  ) : 
    m_Scene(graphicScene),
    m_Title(title),
    m_DrawVertices(drawVertices),
    m_DrawEdges(drawEdges),
    m_DrawFaces(drawFaces),
    m_DrawRays(drawRays),
    m_DrawLines(drawLines),
    m_UseDefaultColor(useDefaultColor),
    m_InverseNormal(inverseNormal),
    m_FlatShading(flatShading)
  {
    initialize();
  } 

  Basic_viewer::~Basic_viewer()
  {
    glDeleteBuffers(NB_GL_BUFFERS, m_VBO);
    glDeleteVertexArrays(NB_VAO_BUFFERS, m_VAO);
    glfwDestroyWindow(m_Window);
    glfwTerminate();
  }

  void Basic_viewer::initialize(bool screenshotOnly)
  {
    m_Window = create_window(m_WindowSize.x(), m_WindowSize.y(), m_Title, screenshotOnly);

    m_AspectRatio = static_cast<float>(m_WindowSize.x())/m_WindowSize.y();

    if (!screenshotOnly)
    {
      // Set event callbacks 
      glfwSetWindowUserPointer(m_Window, this);
      glfwSetKeyCallback(m_Window, key_callback);
      glfwSetCursorPosCallback(m_Window, cursor_callback);
      glfwSetMouseButtonCallback(m_Window, mouse_btn_callback);
      glfwSetScrollCallback(m_Window, scroll_callback);
      glfwSetFramebufferSizeCallback(m_Window, window_size_callback);
    }

    int openglMajorVersion, openglMinorVersion;
    glGetIntegerv(GL_MAJOR_VERSION, &openglMajorVersion);
    glGetIntegerv(GL_MINOR_VERSION, &openglMinorVersion);

    if (openglMajorVersion > 4 || openglMajorVersion == 4 && openglMinorVersion >= 3)
    {
      m_IsOpengl4_3 = true;
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

    check_geometry_feature_availability();

    m_DefaultColorRay = color_to_normalized_vec3(m_Scene.get_default_color_ray());
    m_DefaultColorFace = color_to_normalized_vec3(m_Scene.get_default_color_face());
    m_DefaultColorLine = color_to_normalized_vec3(m_Scene.get_default_color_line());
    m_DefaultColorPoint = color_to_normalized_vec3(m_Scene.get_default_color_point());
    m_DefaultColorSegment = color_to_normalized_vec3(m_Scene.get_default_color_segment());
  }

  void Basic_viewer::show()
  {
    float elapsedTime = 0.0f;
    float lastFrame = 0.0;
    while (!glfwWindowShouldClose(m_Window))
    {
      float currentFrame = static_cast<float>(glfwGetTime());
      m_DeltaTime = currentFrame - lastFrame;
      lastFrame = currentFrame;
      if (m_DeltaTime < 1e-3) m_DeltaTime = 1e-3;

      handle_events(m_DeltaTime);
      if (need_update()) 
      {
        render_scene(m_DeltaTime);
      }
      print_application_state(elapsedTime, m_DeltaTime);
    }
  }

  void Basic_viewer::make_screenshot(const std::string& filePath)
  {
    m_Camera.disable_smoothness();
    draw();
    glfwSwapBuffers(m_Window);
    capture_screenshot(filePath);
  }

  void Basic_viewer::generate_grid(Line_renderer& renderer, float size, int nbSubdivisions) const 
  {
    for (unsigned int i = 0; i <= nbSubdivisions; ++i)
    {
      float pos = float(size * (2.0 * i / nbSubdivisions - 1.0));
      renderer.add_line(vec3f(pos, -size, 0.f), vec3f(pos, size, 0.f));
      renderer.add_line(vec3f(-size, pos, 0.f), vec3f(size, pos, 0.f));
    }
  }

  void Basic_viewer::generate_grid(Line_renderer& renderer, const vec3f& color, float size, int nbSubdivisions) const 
  {
    for (unsigned int i = 0; i <= nbSubdivisions; ++i)
    {
      float pos = float(size * (2.0 * i / nbSubdivisions - 1.0));
      renderer.add_line(vec3f(pos, -size, 0.f), vec3f(pos, size, 0.f), color);
      renderer.add_line(vec3f(-size, pos, 0.f), vec3f(size, pos, 0.f), color);
    }
  }

  void Basic_viewer::compile_shaders()
  {
    const char* PL_VERTEX = m_IsOpengl4_3 ? VERTEX_SOURCE_P_L : VERTEX_SOURCE_P_L_COMP;
    const char* PL_FRAGMENT = m_IsOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    m_ShaderPl = Shader::create(PL_VERTEX, PL_FRAGMENT);

    const char* FACE_VERTEX = m_IsOpengl4_3 ? VERTEX_SOURCE_COLOR : VERTEX_SOURCE_COLOR_COMP;
    const char* FACE_FRAGMENT = m_IsOpengl4_3 ? FRAGMENT_SOURCE_COLOR : FRAGMENT_SOURCE_COLOR_COMP;
    m_ShaderFace = Shader::create(FACE_VERTEX, FACE_FRAGMENT);
    
    const char* PLANE_VERTEX = VERTEX_SOURCE_CLIPPING_PLANE;
    const char* PLANE_FRAGMENT = FRAGMENT_SOURCE_CLIPPING_PLANE;
    m_ShaderPlane = Shader::create(PLANE_VERTEX, PLANE_FRAGMENT);

    const char* SHAPE_VERTEX = VERTEX_SOURCE_SHAPE;
    const char* POINT_GEOMETRY = GEOMETRY_SOURCE_SPHERE;
    m_ShaderSphere = Shader::create(SHAPE_VERTEX, PL_FRAGMENT, POINT_GEOMETRY);

    const char* EDGE_VERTEX = VERTEX_SOURCE_SHAPE;
    const char* EDGE_GEOMETRY = GEOMETRY_SOURCE_CYLINDER;
    const char* EDGE_FRAGMENT = m_IsOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    m_ShaderCylinder = Shader::create(EDGE_VERTEX, PL_FRAGMENT, EDGE_GEOMETRY);

    const char* LINE_VERTEX = VERTEX_SOURCE_LINE_WIDTH;
    const char* LINE_GEOMETRY = GEOMETRY_SOURCE_LINE_WIDTH;
    const char* LINE_FRAGMENT = FRAGMENT_SOURCE_P_L;
    m_ShaderLine = Shader::create(LINE_VERTEX, LINE_FRAGMENT, LINE_GEOMETRY);

    const char* GRID_VERTEX = VERTEX_SOURCE_LINE;
    const char* GRID_GEOMETRY = GEOMETRY_SOURCE_LINE;
    const char* GRID_FRAGMENT = FRAGMENT_SOURCE_LINE;
    m_ShaderGrid = Shader::create(GRID_VERTEX, GRID_FRAGMENT, GRID_GEOMETRY);

    const char* ARROW_VERTEX = VERTEX_SOURCE_LINE;
    const char* ARROW_GEOMETRY = GEOMETRY_SOURCE_ARROW;
    const char* ARROW_FRAGMENT = FRAGMENT_SOURCE_LINE;
    m_ShaderArrow = Shader::create(ARROW_VERTEX, ARROW_FRAGMENT, ARROW_GEOMETRY);

    const char* NORMAL_VERTEX = VERTEX_SOURCE_NORMAL;
    const char* NORMAL_GEOMETRY = GEOMETRY_SOURCE_NORMAL;
    const char* NORMAL_FRAGMENT = m_IsOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    m_ShaderNormal = Shader::create(NORMAL_VERTEX, NORMAL_FRAGMENT, NORMAL_GEOMETRY);

    const char* TRIANGLE_VERTEX = VERTEX_SOURCE_TRIANGLE;
    const char* TRIANGLE_GEOMETRY = GEOMETRY_SOURCE_TRIANGLE;
    const char* TRIANGLE_FRAGMENT = m_IsOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    m_ShaderTriangles = Shader::create(TRIANGLE_VERTEX, TRIANGLE_FRAGMENT, TRIANGLE_GEOMETRY);
  }

  void Basic_viewer::initialize_camera()
  {
    vec3f pmin(
      m_Scene.bounding_box().xmin(),
      m_Scene.bounding_box().ymin(),
      m_Scene.bounding_box().zmin());

    vec3f pmax(
      m_Scene.bounding_box().xmax(),
      m_Scene.bounding_box().ymax(),
      m_Scene.bounding_box().zmax());

    m_BoundingBox = { pmin, pmax };

    m_Camera.lookat(pmin, pmax);  

    if (is_two_dimensional())
    {
      m_Camera.set_constraint_axis(Camera::ConstraintAxis::FORWARD_AXIS);
      m_Camera.set_orthographic();
    }
  }

  void Basic_viewer::initialize_and_load_world_axis()
  {
    // World axis initialization
    m_WorldAxisRenderer.initialize_buffers({ // Use VERTEX_SOURCE_LINE"
      {ShaderDataType::FLOAT3}, // a_Pos"
      {ShaderDataType::FLOAT3}, // a_Color
    });
    m_WorldAxisRenderer.set_width(3.f);
    m_WorldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
    m_WorldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
    m_WorldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitZ(), vec3f(0, 0, 1)); // z-axis
    m_WorldAxisRenderer.load_buffers();

    float cameraSize = m_Camera.get_size() * 0.5;
    // XY grid axis initialization
    m_XYAxisRenderer.initialize_buffers({ // Use VERTEX_SOURCE_LINE
      {ShaderDataType::FLOAT3}, // a_Pos
      {ShaderDataType::FLOAT3}, // a_Color
    });
    m_XYAxisRenderer.set_width(5.f);
    m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
    m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
    m_XYAxisRenderer.load_buffers();

    // XY grid initialization 
    m_XYGridRenderer.initialize_buffers({ // Use VERTEX_SOURCE_LINE
      {ShaderDataType::FLOAT3}, // a_Pos
      {ShaderDataType::FLOAT3}, // a_Color
    });
    m_XYGridRenderer.set_width(2.f);
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitX(), vec3f(.8f, .8f, .8f)); // -x-axis
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitY(), vec3f(.8f, .8f, .8f)); // -y-axis
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitZ(), vec3f(.8f, .8f, .8f)); // -z-axis

    vec3f color(.8f, .8f, .8f);

    generate_grid(m_XYGridRenderer, color, cameraSize);
    
    m_XYGridRenderer.load_buffers();
  }

  void Basic_viewer::load_buffer(int i, int location, const std::vector<float>& vector, int dataCount)
  {
    glBindBuffer(GL_ARRAY_BUFFER, m_VBO[i]);

    glBufferData(GL_ARRAY_BUFFER, vector.size() * sizeof(float), vector.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(location, dataCount, GL_FLOAT, GL_FALSE, dataCount * sizeof(float), nullptr);

    glEnableVertexAttribArray(location);
  }

  void Basic_viewer::load_buffer(int i, int location, int gsEnum, int dataCount)
  {
    const auto& vector = m_Scene.get_array_of_index(gsEnum);
    load_buffer(i, location, vector, dataCount);
  }

  void Basic_viewer::initialize_buffers()
  {
    glGenBuffers(NB_GL_BUFFERS, m_VBO);
    glGenVertexArrays(NB_VAO_BUFFERS, m_VAO);
  }

  void Basic_viewer::load_scene()
  {
    unsigned int bufn = 0;

    std::vector<float> positions, normals, colors; 

    // 1) POINT SHADER

    glBindVertexArray(m_VAO[VAO_POINTS]);
    positions = m_Scene.get_array_of_index(Graphics_scene::POS_POINTS); 
    colors = m_Scene.get_array_of_index(Graphics_scene::COLOR_POINTS);

    load_buffer(bufn++, 0, positions, 3);
    load_buffer(bufn++, 1, colors, 3);

    // 2) SEGMENT SHADER

    glBindVertexArray(m_VAO[VAO_SEGMENTS]);
    positions = m_Scene.get_array_of_index(Graphics_scene::POS_SEGMENTS); 
    colors = m_Scene.get_array_of_index(Graphics_scene::COLOR_SEGMENTS);

    load_buffer(bufn++, 0, positions, 3);
    load_buffer(bufn++, 1, colors, 3);

    // 3) RAYS SHADER

    glBindVertexArray(m_VAO[VAO_RAYS]);
    positions = m_Scene.get_array_of_index(Graphics_scene::POS_RAYS); 
    colors = m_Scene.get_array_of_index(Graphics_scene::COLOR_RAYS);

    load_buffer(bufn++, 0, positions, 3);
    load_buffer(bufn++, 1, colors, 3);

    // 4) LINES SHADER

    glBindVertexArray(m_VAO[VAO_LINES]);
    positions = m_Scene.get_array_of_index(Graphics_scene::POS_LINES); 
    colors = m_Scene.get_array_of_index(Graphics_scene::COLOR_LINES);

    load_buffer(bufn++, 0, positions, 3);
    load_buffer(bufn++, 1, colors, 3);

    // 5) FACE SHADER

    glBindVertexArray(m_VAO[VAO_FACES]);
    positions = m_Scene.get_array_of_index(Graphics_scene::POS_FACES); 
    normals = m_Scene.get_array_of_index(
      m_FlatShading ? Graphics_scene::FLAT_NORMAL_FACES : Graphics_scene::SMOOTH_NORMAL_FACES
    );
    colors = m_Scene.get_array_of_index(Graphics_scene::COLOR_FACES);

    load_buffer(bufn++, 0, positions, 3);
    load_buffer(bufn++, 1, normals, 3);
    load_buffer(bufn++, 2, colors, 3);

    m_AreBuffersInitialized = true;
  }

  CGAL::Plane_3<Basic_viewer::Local_kernel> Basic_viewer::clipping_plane() const
  {
    mat4f CPM = m_ClippingPlane.get_matrix();
    CGAL::Aff_transformation_3<Basic_viewer::Local_kernel> aff(
      CPM(0, 0), CPM(0, 1), CPM(0, 2), CPM(0, 3),
      CPM(1, 0), CPM(1, 1), CPM(1, 2), CPM(1, 3),
      CPM(2, 0), CPM(2, 1), CPM(2, 2), CPM(2, 3)
    );

    CGAL::Plane_3<Local_kernel> p3(0, 0, 1, 0);
    return p3.transform(aff);
  }

  void Basic_viewer::compute_model_view_projection_matrix(const float deltaTime)
  {
    m_Camera.update(deltaTime);
    m_ClippingPlane.update(deltaTime);

    if (m_AnimationController.is_running())
    {
      AnimationKeyFrame animationFrame = m_AnimationController.run();
      m_Camera.set_orientation(animationFrame.orientation);
      m_Camera.set_position(animationFrame.position);
    }

    m_ViewMatrix = m_Camera.view();
    m_ProjectionMatrix = m_Camera.projection(m_WindowSize.x(), m_WindowSize.y());

    m_ViewProjectionMatrix = m_ProjectionMatrix * m_ViewMatrix;
  }

  void Basic_viewer::update_face_uniforms()
  {
    m_ShaderFace->use();

    m_ShaderFace->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderFace->set_mat4f("u_Mv",  m_ViewMatrix.data());
    if (m_UseDefaultColor)
    {
      m_ShaderFace->set_int("u_UseDefaultColor", 1);
      m_ShaderFace->set_vec3f("u_DefaultColor", m_DefaultColorFace.data());
    }
    else 
    {
      m_ShaderFace->set_int("u_UseDefaultColor", 0);
    }

    m_ShaderFace->set_vec4f("u_LightPos",  m_LightPosition.data());
    m_ShaderFace->set_vec4f("u_LightDiff", m_DiffuseColor.data());
    m_ShaderFace->set_vec4f("u_LightSpec", m_SpecularColor.data());
    m_ShaderFace->set_vec4f("u_LightAmb",  m_AmbientColor.data());
    m_ShaderFace->set_float("u_SpecPower", m_Shininess);

    m_ShaderFace->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderFace->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderFace->set_float("u_RenderingTransparency", m_ClippingPlane.get_transparency());
  }

  void Basic_viewer::update_sphere_uniforms()
  {
    m_ShaderSphere->use();

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    if (m_UseDefaultColor)
    {
      m_ShaderSphere->set_int("u_UseDefaultColor", 1);
      m_ShaderSphere->set_vec3f("u_DefaultColor", m_DefaultColorPoint.data());
    }
    else 
    {
      m_ShaderSphere->set_int("u_UseDefaultColor", 0);
    }

    m_ShaderSphere->set_float("u_RenderingMode", static_cast<float>(mode));

    m_ShaderSphere->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderSphere->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderSphere->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderSphere->set_float("u_Radius", m_Camera.get_radius()*m_SizeVertices*0.001);
  }

  void Basic_viewer::update_cylinder_uniforms()
  {
    m_ShaderCylinder->use();

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    if (m_UseDefaultColor)
    {
      m_ShaderCylinder->set_int("u_UseDefaultColor", 1);
      m_ShaderCylinder->set_vec3f("u_DefaultColor", m_DefaultColorSegment.data());
    }
    else 
    {
      m_ShaderCylinder->set_int("u_UseDefaultColor", 0);
    }

    m_ShaderCylinder->set_float("u_RenderingMode", static_cast<float>(mode));

    m_ShaderCylinder->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderCylinder->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderCylinder->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderCylinder->set_float("u_Radius", m_Camera.get_radius()*m_SizeEdges*0.001);
  }

  void Basic_viewer::update_pl_uniforms(const vec3f& defaultColor)
  {
    m_ShaderPl->use();

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    m_ShaderPl->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderPl->set_float("u_PointSize", m_Camera.get_radius()*m_SizeVertices*0.1);
    m_ShaderPl->set_int("u_IsOrthographic", static_cast<int>(m_Camera.is_orthographic()));
    if (m_UseDefaultColor)
    {
      m_ShaderPl->set_int("u_UseDefaultColor", 1);
      m_ShaderPl->set_vec3f("u_DefaultColor", defaultColor.data());
    }
    else 
    {
      m_ShaderPl->set_int("u_UseDefaultColor", 0);
    }

    m_ShaderPl->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderPl->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderPl->set_float("u_RenderingMode", static_cast<float>(mode));
  }

  void Basic_viewer::update_line_uniforms(float size, const vec3f& defaultColor)
  {
    m_ShaderLine->use(); 

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    float viewport[2] = { m_WindowSize.x(), m_WindowSize.y() }; 

    m_ShaderLine->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderLine->set_float("u_PointSize", size);
    m_ShaderLine->set_int("u_IsOrthographic", static_cast<int>(m_Camera.is_orthographic()));
    if (m_UseDefaultColor)
    {
      m_ShaderLine->set_int("u_UseDefaultColor", 1);
      m_ShaderLine->set_vec3f("u_DefaultColor", defaultColor.data());
    }
    else 
    {
      m_ShaderLine->set_int("u_UseDefaultColor", 0);
    }
    
    m_ShaderLine->set_vec2f("u_Viewport", &viewport[0]);

    m_ShaderLine->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderLine->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderLine->set_float("u_RenderingMode", static_cast<float>(mode));
  }

  void Basic_viewer::update_clipping_uniforms()
  {
    mat4f clippingModelMatrix = m_ClippingPlane.get_matrix();

    m_PointPlane = clippingModelMatrix * vec4f(0, 0, 0, 1);
    m_ClipPlane = clippingModelMatrix * vec4f(0, 0, 1, 0);

    m_ShaderPlane->use();
    m_ShaderPlane->set_mat4f("u_Vp", m_ViewProjectionMatrix.data());
    m_ShaderPlane->set_mat4f("u_M",  clippingModelMatrix.data());
  }

  void Basic_viewer::update_world_axis_uniforms()
  {
    mat4f view = m_ViewMatrix;

    // we only want the rotation part of the view matrix  
    mat3f rotation = view.block<3,3>(0,0);

    mat4f rotation4x4 = mat4f::Identity();
    rotation4x4.block<3,3>(0,0) = rotation;

    float halfWidth = m_AspectRatio * 0.1f;
    float halfHeight = 0.1f;
    mat4f projection = utils::ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0f, 1.0f);

    mat4f translation = transform::translation(vec3f(halfWidth - 0.1f*m_AspectRatio, halfHeight - 0.1f, 0.0f));

    mat4f mvp = projection * rotation4x4 * translation;

    m_ShaderArrow->use();
    m_ShaderArrow->set_mat4f("u_Mvp", mvp.data()); 
    m_ShaderArrow->set_float("u_SceneRadius", 1.0f); 
  }

  void Basic_viewer::update_XY_axis_uniforms()
  {
    m_ShaderArrow->use();
    m_ShaderArrow->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
    m_ShaderArrow->set_float("u_SceneRadius", m_Camera.get_radius()); 
  }

  void Basic_viewer::update_XY_grid_uniforms()
  {
    m_ShaderGrid->use();
    m_ShaderGrid->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());
  }

  void Basic_viewer::update_normals_uniforms()
  {
    m_ShaderNormal->use();

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    vec4f color = color_to_normalized_vec4(m_DefaultColorNormal);
    m_ShaderNormal->set_mat4f("u_Mv", m_ViewMatrix.data());
    if (m_UseDefaultColorNormal)
    {
      m_ShaderNormal->set_int("u_UseDefaultColor", 1);
      m_ShaderNormal->set_vec3f("u_DefaultColor", color.data());
    }
    else
    {
      m_ShaderNormal->set_int("u_UseDefaultColor", 0);
    }
    m_ShaderNormal->set_mat4f("u_Projection", m_ProjectionMatrix.data());
    m_ShaderNormal->set_float("u_Factor", m_NormalHeightFactor);
    m_ShaderNormal->set_float("u_SceneRadius", m_Camera.get_radius());
    if (m_DisplayFaceNormal)
    {
      m_ShaderNormal->set_int("u_DisplayFaceNormal", 1);
    }
    else
    {
      m_ShaderNormal->set_int("u_DisplayFaceNormal", 0);
    }

    m_ShaderNormal->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderNormal->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderNormal->set_float("u_RenderingMode", static_cast<float>(mode));
  }

  void Basic_viewer::update_triangles_uniforms()
  {
    m_ShaderTriangles->use();

    bool half = m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    m_ShaderTriangles->set_mat4f("u_Mvp", m_ViewProjectionMatrix.data());

    m_ShaderTriangles->set_vec4f("u_ClipPlane",  m_ClipPlane.data());
    m_ShaderTriangles->set_vec4f("u_PointPlane", m_PointPlane.data());
    m_ShaderTriangles->set_float("u_RenderingMode", static_cast<float>(mode));
  }

  bool Basic_viewer::need_update() const
  {
    return m_Camera.need_update() 
        || m_ClippingPlane.need_update() 
        || m_AnimationController.is_running()
        || has_active_actions() 
        ; 
  }

  void Basic_viewer::render_scene(const float deltaTime)
  {
    draw(deltaTime);
    glfwSwapBuffers(m_Window);
  }

  void Basic_viewer::draw(const float deltaTime)
  {
    if (!m_AreBuffersInitialized)
    {
      load_scene();
    }

    glClearColor(1.0f, 1.0f, 1.0f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);

    compute_model_view_projection_matrix(deltaTime);

    if (m_DrawRays)
    {
      draw_rays();
    }

    if (m_DrawLines)
    {
      draw_lines();
    }
   
    if (m_DrawEdges && !m_DrawMeshTriangles)
    {
      draw_edges();
    }
    
    if (m_DrawVertices)
    {
      draw_vertices();
    }

    if (clipping_plane_enable())
    {
      render_clipping_plane();
    }

    if (m_DrawNormals)
    {
      draw_normals();
    }  
        
    if (m_DrawFaces)
    {
      draw_faces();
    }

    if (m_DrawMeshTriangles)
    {
      draw_triangles();
    }
    
    if (m_DrawWorldAxis)  
    {
      draw_world_axis();
    }

    if (m_DrawXYGrid) 
    { 
      draw_xy_grid();
    }
  }

  inline
  vec4f Basic_viewer::color_to_normalized_vec4(const CGAL::IO::Color& c) const
  {
    return { static_cast<float>(c.red()) / 255, static_cast<float>(c.green()) / 255, static_cast<float>(c.blue()) / 255, 1.0f };
  }

  inline
  vec3f Basic_viewer::color_to_normalized_vec3(const CGAL::IO::Color& c) const
  {
    return { static_cast<float>(c.red()) / 255, static_cast<float>(c.green()) / 255, static_cast<float>(c.blue()) / 255 };
  }

  void Basic_viewer::draw_faces()
  {
    update_face_uniforms();

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2.0, 2.0); 
    glDepthFunc(GL_LESS);
    if (m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF)
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
    } 
    else // Not CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF
    {
      if (m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
          m_DisplayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        // 1. draw solid HALF
        draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);
      } 
      else
      {
        // 1. draw solid FOR ALL
        draw_faces_bis(RenderingMode::DRAW_ALL);
      } 
    }
    glDisable(GL_POLYGON_OFFSET_FILL); 
  }

  void Basic_viewer::draw_faces_bis(RenderingMode mode)
  {
    m_ShaderFace->set_float("u_RenderingMode", static_cast<float>(mode));

    glBindVertexArray(m_VAO[VAO_FACES]);
    glDrawArrays(GL_TRIANGLES, 0, m_Scene.number_of_elements(Graphics_scene::POS_FACES));
  }

  void Basic_viewer::draw_rays()
  {
    update_pl_uniforms(m_DefaultColorRay);

    m_ShaderPl->set_float("u_RenderingMode", static_cast<float>(RenderingMode::DRAW_ALL));

    glLineWidth(m_SizeRays);
    glBindVertexArray(m_VAO[VAO_RAYS]);
    glDrawArrays(GL_LINES, 0, m_Scene.number_of_elements(Graphics_scene::POS_RAYS));
    glLineWidth(1.0);
  }

  void Basic_viewer::draw_vertices()
  {
    update_pl_uniforms(m_DefaultColorPoint);
    if (m_DrawSphereVertex && m_GeometryFeatureEnabled)
    {
      update_sphere_uniforms();
    }

    glBindVertexArray(m_VAO[VAO_POINTS]);
    glDrawArrays(GL_POINTS, 0, m_Scene.number_of_elements(Graphics_scene::POS_POINTS));
  }

  void Basic_viewer::draw_lines()
  {
    update_pl_uniforms(m_DefaultColorLine);

    m_ShaderPl->set_float("u_RenderingMode", static_cast<float>(RenderingMode::DRAW_ALL));

    glLineWidth(m_SizeLines);
    glBindVertexArray(m_VAO[VAO_LINES]);
    glDrawArrays(GL_LINES, 0, m_Scene.number_of_elements(Graphics_scene::POS_LINES));
    glLineWidth(1.0);
  }

  void Basic_viewer::draw_edges()
  {
    update_line_uniforms(m_SizeEdges, m_DefaultColorSegment);
    if (m_DrawCylinderEdge && m_GeometryFeatureEnabled)
    {
      update_cylinder_uniforms();
    }

    glDepthFunc(GL_LEQUAL);
    glBindVertexArray(m_VAO[VAO_SEGMENTS]);
    glDrawArrays(GL_LINES, 0, m_Scene.number_of_elements(Graphics_scene::POS_SEGMENTS));
  }


  void Basic_viewer::draw_world_axis() 
  {
    update_world_axis_uniforms();
    glDepthFunc(GL_LEQUAL);

    int &w = m_WindowSize.x();
    int &h = m_WindowSize.y();

    glViewport(w - w / 5, h - h / 5, w / 5, h / 5);
    m_WorldAxisRenderer.draw();

    // Restore the main viewport
    glViewport(0, 0, w, h);
  }

  void Basic_viewer::draw_xy_grid() 
  { 
    glDepthFunc(GL_LEQUAL);

    update_XY_grid_uniforms();
    m_XYGridRenderer.draw();
    update_XY_axis_uniforms();
    m_XYAxisRenderer.draw();
  }

  void Basic_viewer::draw_normals()
  {
    update_normals_uniforms();

    glDepthFunc(GL_LEQUAL);
    glBindVertexArray(m_VAO[VAO_FACES]);
    glLineWidth(m_SizeNormals);
    glDrawArrays(GL_TRIANGLES, 0, m_Scene.number_of_elements(Graphics_scene::POS_FACES));
  }

  void Basic_viewer::draw_triangles()
  {
    update_triangles_uniforms();

    glDepthFunc(GL_LEQUAL);
    glBindVertexArray(m_VAO[VAO_FACES]);
    glDrawArrays(GL_TRIANGLES, 0, m_Scene.number_of_elements(Graphics_scene::POS_FACES));
  }

  void Basic_viewer::initialize_and_load_clipping_plane()
  {
    float size = ((m_Scene.bounding_box().xmax() - m_Scene.bounding_box().xmin()) +
                  (m_Scene.bounding_box().ymax() - m_Scene.bounding_box().ymin()) +
                  (m_Scene.bounding_box().zmax() - m_Scene.bounding_box().zmin()));

    const unsigned int NB_SUBDIVISIONS = 30;

    m_ClippingPlane.initialize_buffers({ // Use VERTEX_SOURCE_CLIPPING_PLANE
      {ShaderDataType::FLOAT3}, // a_Pos
    });
    m_ClippingPlane.set_width(0.2f);
    m_ClippingPlane.add_line({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
    generate_grid(m_ClippingPlane, size, NB_SUBDIVISIONS);
    m_ClippingPlane.load_buffers();

    m_ClippingPlane.set_size(m_Camera.get_size());  
  }

  void Basic_viewer::render_clipping_plane()
  {
    if (!m_IsOpengl4_3)
    {
      return;
    }

    update_clipping_uniforms();
    if (m_DrawClippingPlane)
    {
      m_ClippingPlane.draw();
    }
  }

  void Basic_viewer::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_key_event(key, scancode, action, mods);
  }

  void Basic_viewer::cursor_callback(GLFWwindow* window, double xpos, double ypo)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    viewer->on_cursor_event(xpos, ypo, windowWidth, windowHeight);
  }

  void Basic_viewer::mouse_btn_callback(GLFWwindow* window, int button, int action, int mods)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_mouse_button_event(button, action, mods);
  }

  void Basic_viewer::window_size_callback(GLFWwindow* window, int width, int height)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->m_WindowSize = {width, height};

    viewer->m_AspectRatio = static_cast<float>(width) / height;
    glViewport(0, 0, width, height);

    viewer->render_scene(viewer->m_DeltaTime);
  }

  void Basic_viewer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_scroll_event(xoffset, yoffset);
  }
 
  void Basic_viewer::print_application_state(float& elapsedTime, const float deltaTime)
  {
    elapsedTime += deltaTime;
    if (elapsedTime * 1000 > 100) // update terminal display each 100ms
    {
      elapsedTime = 0.0f;
      if (m_PrintApplicationState)
      {
        std::cout << "\33[2K"  
                  << std::round(1 / deltaTime)      << " fps    " 
                  << deltaTime*1000                 << " ms\n\33[2K" 
                  << "Camera type: "                << (m_Camera.is_orbiter() ? "ORBITER" : "FREE_FLY")               << "    " 
                  << "Camera mode: "                << (m_Camera.is_orthographic() ? "ORTHOGRAPHIC" : "PERSPECTIVE")  << "    " 
                  << "FOV: "                        << m_Camera.get_fov()                                             << " \n\33[2K" 
                  << "Camera translation speed: "   << m_Camera.get_translation_speed()                        << "    " 
                  << "Camera rotation speed: "      << std::round(m_Camera.get_rotation_speed())               << "    "
                  << "Camera constraint axis: "     << m_Camera.get_constraint_axis_str()                      << "\n\33[2K"     
                  << "CP translation speed: "       << m_ClippingPlane.get_translation_speed()                 << "    "            
                  << "CP rotation speed: "          << std::round(m_ClippingPlane.get_rotation_speed())        << "    "     
                  << "CP constraint axis: "         << m_ClippingPlane.get_constraint_axis_str()               << "\n\33[2K" 
                  << "Reversed normals : "          << (m_InverseNormal ? "TRUE" : "FALSE")                    << "    " 
                  << "Face normal : "               << (m_DisplayFaceNormal ? "TRUE" : "FALSE")                << "    " 
                  << "Shading : "                   << (m_FlatShading ? "FLAT" : "SMOOTH")                          << "\n\33[2K"
                  << "Draw faces: "                 << (m_DrawFaces ? "TRUE" : "FALSE")                        << "    " 
                  << "Draw edges: "                 << (m_DrawEdges ? "TRUE" : "FALSE")                        << "    " 
                  << "Draw vertices: "              << (m_DrawVertices ? "TRUE" : "FALSE")                     << "\n\33[2K"
                  << "Draw lines: "                 << (m_DrawLines ? "TRUE" : "FALSE")                        << "    " 
                  << "Draw rays: "                  << (m_DrawRays ? "TRUE" : "FALSE")                         << "    " 
                  << "Draw normals: "               << (m_DrawNormals ? "TRUE" : "FALSE")                      << "\n\33[2K"  
                  << "Draw edges as cylinder: "     << (m_DrawCylinderEdge ? "TRUE" : "FALSE")                 << "    " 
                  << "Draw vertices as sphere: "    << (m_DrawSphereVertex ? "TRUE" : "FALSE")                 << "\n\33[2K"    
                  << "Use default color: "          << (m_UseDefaultColor ? "TRUE" : "FALSE")                  << "    " 
                  << "Number of key frames: "       << m_AnimationController.number_of_key_frames()            << "\n\33[2K"     
                  << "Light color: ("               << m_AmbientColor.x()                                      << ", " 
                                                    << m_AmbientColor.y()                                      << ", " 
                                                    << m_AmbientColor.z()                                      << ")\n\33[2K"     
                  << "Size of vertices: "           << m_SizeVertices << "    Size of edges: " <<  m_SizeEdges << "    "            
                  << "\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\r"   << std::flush;
      }
    }
  }

  void Basic_viewer::start_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case ROTATE_CLIPPING_PLANE:
    case TRANSLATE_CLIPPING_PLANE:
    case TRANSLATE_CP_ALONG_CAMERA_DIRECTION:
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      // glfwSetInputMode(m_Window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      break;
    }
  }
 
  void Basic_viewer::exit_app() 
  {
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\nAPPLICATION IS CLOSING" << std::endl;
    exit(EXIT_SUCCESS);
  }

  void Basic_viewer::check_geometry_feature_availability()
  {
    int maxGeometryOutputVertices = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES, &maxGeometryOutputVertices);
    int maxGeometryOutputComponents = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS, &maxGeometryOutputComponents);

    if (maxGeometryOutputVertices < 128 || maxGeometryOutputComponents < 1024)
    {
      std::cout << "Cylinder edge and sphere vertex feature disabled! (maxGeometryOutputVertices=" << maxGeometryOutputVertices << ", maxGeometryOutputComponents=" << maxGeometryOutputComponents << ")\n";
      m_GeometryFeatureEnabled = false; 
    }
  }

  void Basic_viewer::double_click_event(int btn)
  {
    if (m_Camera.is_orbiter()) 
    {
      if (btn == GLFW_MOUSE_BUTTON_RIGHT)
      {
        if (is_key_pressed(m_Window, GLFW_KEY_LEFT_CONTROL)) 
        {
          m_Camera.reset_orientation();
        }
        m_Camera.reset_position();
      }
      else if (btn == GLFW_MOUSE_BUTTON_LEFT)
      {
        if (is_key_pressed(m_Window, GLFW_KEY_LEFT_CONTROL))
        {
          m_ClippingPlane.align_to_direction(m_Camera.get_forward());
        }
        else if (is_key_pressed(m_Window, GLFW_KEY_LEFT_SHIFT))
        {
          m_Camera.align_to_plane(m_ClippingPlane.get_normal());
        }
        else
        {
          m_Camera.align_to_nearest_axis();
        }
      }
      else if (btn == GLFW_MOUSE_BUTTON_MIDDLE)
      {
        m_Camera.reset_size();
      }
    }
  }

  void Basic_viewer::scroll_event(const float deltaTime)
  {
    float yoffset = get_scroll_yOffset() / m_AspectRatio;

    if (is_key_pressed(m_Window, GLFW_KEY_LEFT_SHIFT))
    {
      m_Camera.increase_zoom_smoothness(yoffset);
    }
    else if (is_key_pressed(m_Window, GLFW_KEY_Z) && !m_Camera.is_orthographic())
    {
      m_Camera.increase_fov(yoffset);
      m_ClippingPlane.set_size(m_Camera.get_size());
    }
    else if (is_key_pressed(m_Window, GLFW_KEY_LEFT_CONTROL) && m_DisplayMode != DisplayMode::CLIPPING_PLANE_OFF)
    {
        m_ClippingPlane.translation(8.f * yoffset * deltaTime);
    }
    else 
    {
      m_Camera.move(8.f * yoffset * deltaTime);
      m_ClippingPlane.set_size(m_Camera.get_size());
    }
  }
 
  void Basic_viewer::change_pivot_point() 
  {
    auto [mouseX, mouseY] = get_mouse_position();

    vec2f nc = utils::normalized_coordinates({mouseX, mouseY}, m_WindowSize.x(), m_WindowSize.y());

    vec3f cameraPosition = m_Camera.get_position(); 
    float cameraX = cameraPosition.x();
    float cameraY = cameraPosition.y();

    float xValue = nc.x() + cameraX; 
    float yValue = nc.y() + cameraY;

    // vec4f test = {cameraPosition.x(), cameraPosition.y(), cameraPosition.z(), 1.0};

    // vec4f vppos = transform::viewport(m_WindowSize.x(), m_WindowSize.y) * 

    if (utils::inside_bounding_box_2d({xValue, yValue}, {m_BoundingBox.first.x(), m_BoundingBox.first.y()}, {m_BoundingBox.second.x(), m_BoundingBox.second.y()})) 
    {
      std::cout << "INSIDE\n";
      m_Camera.set_center({nc.x(), nc.y(), 0});
    }
    else 
    {
      m_Camera.set_center(utils::center(m_BoundingBox.first, m_BoundingBox.second));
      std::cout << "OUTSIDE\n";
    }

    // std::cout << "Mouse position : " << nc.x() << " " << nc.y() << ", Camera position : " << cameraPosition.transpose() << "\n";
    // std::cout << "BBOX : " << m_BoundingBox.first.transpose() << " " << m_BoundingBox.second.transpose() << "\n";
  }

  void Basic_viewer::action_event(int action, const float deltaTime)
  {
    switch (action)
    {
    /*APPLICATION*/
    case EXIT: 
      exit_app();
    case PRINT_APPLICATION_STATE:
      m_PrintApplicationState = !m_PrintApplicationState; 
      std::cout << "\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K"
                << "\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\r" << std::flush;
      break;
    /*WINDOW*/
    case FULLSCREEN:
      fullscreen();
      break;
    case SCREENSHOT:
      capture_screenshot("./screenshot.png");
      break;
    /*SCENE*/
    case NORMALS_DISPLAY:
      m_DrawNormals = !m_DrawNormals;
      break;
    case FACE_NORMALS_DISPLAY:
      m_DisplayFaceNormal = !m_DisplayFaceNormal; 
      break;
    case TRIANGLES_DISPLAY: 
      m_DrawMeshTriangles = !m_DrawMeshTriangles;
      break;
    case VERTICES_DISPLAY:
      m_DrawVertices = !m_DrawVertices;
      break;
    case SPHERE_VERTEX_DISPLAY:
      m_DrawSphereVertex = !m_DrawSphereVertex;
      break;
    case FACES_DISPLAY:
      m_DrawFaces = !m_DrawFaces;
      break;
    case EDGES_DISPLAY:
      m_DrawEdges = !m_DrawEdges;
      break;
    case CYLINDER_EDGE_DISPLAY:
      m_DrawCylinderEdge = !m_DrawCylinderEdge;
      break;
    case SHADING_MODE:
      m_FlatShading = !m_FlatShading;
      m_AreBuffersInitialized = false;
      break;
    case INVERSE_NORMAL:
      m_InverseNormal = !m_InverseNormal;
      m_Scene.reverse_all_normals();
      m_AreBuffersInitialized = false;
      break;
    case MONO_COLOR:
      m_UseDefaultColor = !m_UseDefaultColor;
      break;
    case NORMALS_MONO_COLOR: 
      m_UseDefaultColorNormal = !m_UseDefaultColorNormal;
      break;
    case INC_EDGES_SIZE:
      increase_size_edge(deltaTime);
      break;
    case DEC_EDGES_SIZE:
      decrease_size_edge(deltaTime);
      break;
    case INC_POINTS_SIZE:
      increase_size_vertex(deltaTime);
      break;
    case DEC_POINTS_SIZE:
      decrease_size_vertex(deltaTime);
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
    case WORLD_AXIS_DISPLAY:
      m_DrawWorldAxis = !m_DrawWorldAxis;
      break;
    case XY_GRID_DISPLAY:
      m_DrawXYGrid = !m_DrawXYGrid;
      break;
    /*CAMERA*/
    case UP:
      m_Camera.move_up(deltaTime);
      break;
    case DOWN:
      m_Camera.move_down(deltaTime);
      break;
    case LEFT:
      m_Camera.move_left(deltaTime);
      break;
    case RIGHT:
      m_Camera.move_right(deltaTime);
      break;
    case FORWARD:
      m_Camera.move(deltaTime);
      break;
    case BACKWARDS:
      m_Camera.move(-deltaTime);
      break;
    case INC_CAMERA_ROTATION_SMOOTHNESS:
      m_Camera.increase_rotation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SMOOTHNESS:
      m_Camera.decrease_rotation_smoothness(deltaTime);
      break;
    case INC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_Camera.increase_translation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_Camera.decrease_translation_smoothness(deltaTime);
      break;
    case SWITCH_CAMERA_MODE:
      m_Camera.toggle_mode();
      break;
    case SWITCH_CAMERA_TYPE:
      m_Camera.toggle_type();
      break;
    case SWITCH_CAMERA_CONSTRAINT_AXIS: 
      m_Camera.switch_constraint_axis();
      break;
    case INC_CAMERA_TRANSLATION_SPEED:
      m_Camera.increase_translation_speed(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SPEED:
      m_Camera.decrease_translation_speed(deltaTime);
      break;
    case INC_CAMERA_ROTATION_SPEED:
      m_Camera.increase_rotation_speed(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SPEED:
      m_Camera.decrease_rotation_speed(deltaTime);
      break;
    case ROTATE_CAMERA:
      rotate_camera();
      break;
    case TRANSLATE_CAMERA:
      translate_camera(deltaTime);
      break;
    case RESET_CAMERA_AND_CP:
      reset_camera_and_clipping_plane();
      break;
    case CHANGE_PIVOT_POINT:
      change_pivot_point();
      break;
    /*CLIPPING PLANE*/
    case INC_CP_TRANSLATION_SPEED:
      m_ClippingPlane.increase_translation_speed(deltaTime);
      break;
    case DEC_CP_TRANSLATION_SPEED:
      m_ClippingPlane.decrease_translation_speed(deltaTime);
      break;
    case INC_CP_ROTATION_SPEED:
      m_ClippingPlane.increase_rotation_speed(deltaTime);
      break;
    case DEC_CP_ROTATION_SPEED:
      m_ClippingPlane.decrease_rotation_speed(deltaTime);
      break;
    case SWITCH_CLIPPING_PLANE_DISPLAY:
      m_DrawClippingPlane = !m_DrawClippingPlane;
      break;
    case SWITCH_CLIPPING_PLANE_MODE:
      switch_display_mode();
      break;
    case ROTATE_CLIPPING_PLANE:
      rotate_clipping_plane();
      break;
    case TRANSLATE_CLIPPING_PLANE:
      translate_clipping_plane(deltaTime);
      break;
    case TRANSLATE_CP_ALONG_CAMERA_DIRECTION:
      translate_clipping_plane(deltaTime, true);
      break;
    case SWITCH_CP_CONSTRAINT_AXIS:
      m_ClippingPlane.switch_constraint_axis();
      break;
    case RESET_CLIPPING_PLANE: 
      reset_clipping_plane();
      break;
    /*ANIMATION*/
    case SAVE_KEY_FRAME:
      save_key_frame();
      break;
    case RUN_OR_STOP_ANIMATION:
      run_or_stop_animation();
      break;
    case CLEAR_ANIMATION: 
      m_AnimationController.clear_buffer();
      break;
    }
  }
 
  void Basic_viewer::increase_size_edge(const float deltaTime)
  {
    float upperBound = 20.0;  
    m_SizeEdges = std::min(upperBound, m_SizeEdges + 10.0f*deltaTime);
  }

  void Basic_viewer::decrease_size_edge(const float deltaTime)
  {
    float lowerBound = 0.01;
    m_SizeEdges = std::max(lowerBound, m_SizeEdges - 10.0f*deltaTime);
  }

  void Basic_viewer::increase_size_vertex(const float deltaTime)
  {
    float upperBound = 50.0;  
    m_SizeVertices = std::min(upperBound, m_SizeVertices + 10.0f*deltaTime);
  }

  void Basic_viewer::decrease_size_vertex(const float deltaTime)
  {
    float lowerBound = 0.01;
    m_SizeVertices = std::max(lowerBound, m_SizeVertices - 10.0f*deltaTime);
  }

  void Basic_viewer::increase_red_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_AmbientColor.x() += speed;
    if (m_AmbientColor.x() > 1.f)
      m_AmbientColor.x() = 1.f;
    if (m_AmbientColor.x() < 0.f)
      m_AmbientColor.x() = 0.f;

    m_DiffuseColor.x() += speed;
    if (m_DiffuseColor.x() > 1.f)
      m_DiffuseColor.x() = 1.f;
    if (m_DiffuseColor.x() < 0.f)
      m_DiffuseColor.x() = 0.f;

    m_SpecularColor.x() += speed;
    if (m_SpecularColor.x() > 1.f)
      m_SpecularColor.x() = 1.f;
    if (m_SpecularColor.x() < 0.f)
      m_SpecularColor.x() = 0.f;
  }
 
  void Basic_viewer::increase_green_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_AmbientColor.y() += speed;
    if (m_AmbientColor.y() > 1.f)
      m_AmbientColor.y() = 1.f;
    if (m_AmbientColor.y() < 0.f)
      m_AmbientColor.y() = 0.f;

    m_DiffuseColor.y() += speed;
    if (m_DiffuseColor.y() > 1.f)
      m_DiffuseColor.y() = 1.f;
    if (m_DiffuseColor.y() < 0.f)
      m_DiffuseColor.y() = 0.f;

    m_SpecularColor.y() += speed;
    if (m_SpecularColor.y() > 1.f)
      m_SpecularColor.y() = 1.f;
    if (m_SpecularColor.y() < 0.f)
      m_SpecularColor.y() = 0.f;
  }
 
  void Basic_viewer::increase_blue_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_AmbientColor.z() += speed;
    if (m_AmbientColor.z() > 1.f)
      m_AmbientColor.z() = 1.f;
    if (m_AmbientColor.z() < 0.f)
      m_AmbientColor.z() = 0.f;

    m_DiffuseColor.z() += speed;
    if (m_DiffuseColor.z() > 1.f)
      m_DiffuseColor.z() = 1.f;
    if (m_DiffuseColor.z() < 0.f)
      m_DiffuseColor.z() = 0.f;

    m_SpecularColor.z() += speed;
    if (m_SpecularColor.z() > 1.f)
      m_SpecularColor.z() = 1.f;
    if (m_SpecularColor.z() < 0.f)
      m_SpecularColor.z() = 0.f;
  }

 
  void Basic_viewer::increase_light_all(const float deltaTime)
  {
    increase_red_component(deltaTime);
    increase_green_component(deltaTime);
    increase_blue_component(deltaTime);
  }
 
  void Basic_viewer::save_key_frame() 
  {
    if (m_Camera.is_orbiter())
    {
      m_AnimationController.add_key_frame(
        m_Camera.get_position(),
        m_Camera.get_orientation()
      );
    }
  }
 
  void Basic_viewer::run_or_stop_animation() 
  {
    if (m_Camera.is_orbiter())
    {
      if (m_AnimationController.is_running()) 
      {
        m_AnimationController.stop(m_AnimationController.get_frame());
      }
      else 
      {
        m_AnimationController.start();
      }
    }
  }

  void Basic_viewer::reset_camera_and_clipping_plane()
  {
    m_Camera.reset_all();
    m_ClippingPlane.reset_all();
    m_ClippingPlane.set_size(m_Camera.get_size());
  }

  void Basic_viewer::reset_clipping_plane()
  {
    m_ClippingPlane.reset_all();
    m_ClippingPlane.set_size(m_Camera.get_size());
  }

  void Basic_viewer::switch_display_mode()
  {
    switch(m_DisplayMode) 
    {
    case DisplayMode::CLIPPING_PLANE_OFF:
      m_DisplayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF:
      m_DisplayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF: 
      m_DisplayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY: 
      m_DisplayMode = DisplayMode::CLIPPING_PLANE_OFF; 
      break;
    }
  }

  void Basic_viewer::rotate_clipping_plane()
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    m_ClippingPlane.set_right_axis(m_Camera.get_right());
    m_ClippingPlane.set_up_axis(m_Camera.get_up());

    m_ClippingPlane.rotation(
      mouseDeltaX / m_AspectRatio, 
      mouseDeltaY / m_AspectRatio
    );
  }

  void Basic_viewer::translate_clipping_plane(const float deltaTime, bool useCameraForward)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    float deltaX = deltaTime * mouseDeltaX / m_AspectRatio;
    float deltaY = deltaTime * mouseDeltaY / m_AspectRatio;

    if (useCameraForward)
    {
      vec3f forwardDirection = m_Camera.get_forward();

      float s = abs(deltaY) > abs(deltaX) ? -deltaY : deltaX;

      m_ClippingPlane.translation(forwardDirection, s);
    }
    else 
    {
      m_ClippingPlane.translation(-deltaX, deltaY);
    }
  }

  void Basic_viewer::rotate_camera()
  {
    auto mouseDelta = get_mouse_delta();

    if (m_Camera.get_constraint_axis_str() == "Forward")
    {
      mouseDelta = get_roll_mouse_delta();
    }

    float mouseDeltaX = mouseDelta.first;
    float mouseDeltaY = mouseDelta.second; 
    
    m_Camera.rotation(
      mouseDeltaX / m_AspectRatio, 
      mouseDeltaY / m_AspectRatio
    );
  }

  void Basic_viewer::translate_camera(const float deltaTime)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    m_Camera.translation(
      deltaTime * -mouseDeltaX / m_AspectRatio,
      deltaTime * mouseDeltaY / m_AspectRatio
    );
  }

  void Basic_viewer::fullscreen()
  {
    m_IsFullscreen = !m_IsFullscreen;

    int count;
    GLFWmonitor *monitor = glfwGetMonitors(&count)[0];
    const GLFWvidmode *mode = glfwGetVideoMode(monitor);

    if (m_IsFullscreen) 
    {
      m_OldWindowSize = m_WindowSize;

#if !defined(GLFW_USE_WAYLAND)
      glfwGetWindowPos(m_Window, &m_OldWindowPosition.x(), &m_OldWindowPosition.y()); 
#endif 
      glfwSetWindowMonitor(m_Window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
      glViewport(0, 0, mode->width, mode->height);

      m_WindowSize.x() = mode->width;
      m_WindowSize.y() = mode->height;
    }
    else 
    {
      m_WindowSize = m_OldWindowSize;
      glfwSetWindowMonitor(m_Window, nullptr, m_OldWindowPosition.x(), m_OldWindowPosition.y(), m_WindowSize.x(), m_WindowSize.y(), mode->refreshRate);
      glViewport(0, 0, m_WindowSize.x(), m_WindowSize.y());
    }
  }

  void Basic_viewer::capture_screenshot(const std::string &filepath)
  {
    // https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/ (thanks)
    // https://github.com/nothings/stb/
    // The stb lib used here is from glfw/deps

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);

    const GLsizei NB_CHANNELS = 4;
    GLsizei stride = NB_CHANNELS * m_WindowSize.x();
    stride += (stride % 4) ? (4 - stride % 4) : 0; // stride must be a multiple of 4
    GLsizei bufferSize = stride * m_WindowSize.y();

    std::vector<char> buffer(bufferSize);
    m_ShaderFace->use(); 
    glPixelStorei(GL_PACK_ALIGNMENT, 4);

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, m_WindowSize.x(), m_WindowSize.y(), GL_RGBA, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath.data(), m_WindowSize.x(), m_WindowSize.y(), NB_CHANNELS, buffer.data(), stride);
  }
 
  void Basic_viewer::end_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      glfwSetInputMode(m_Window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
      break;
    }
  }

  void Basic_viewer::initialize_keys_actions()
  {
    /*APPLICATION*/
    add_keyboard_action({GLFW_KEY_ESCAPE                  }, InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_Q, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_T                       }, InputMode::RELEASE, PRINT_APPLICATION_STATE);

    /*WINDOW*/
    add_keyboard_action({GLFW_KEY_ENTER, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, FULLSCREEN);
    add_keyboard_action({GLFW_KEY_F2                      }, InputMode::RELEASE, SCREENSHOT);

    /*SCENE*/

    add_keyboard_action({GLFW_KEY_N, GLFW_KEY_LEFT_SHIFT}, InputMode::RELEASE, FACE_NORMALS_DISPLAY);
    add_keyboard_action({GLFW_KEY_N, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, NORMALS_DISPLAY);
    add_keyboard_action({GLFW_KEY_M, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, NORMALS_MONO_COLOR);

    add_keyboard_action({GLFW_KEY_T, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, TRIANGLES_DISPLAY);

    add_keyboard_action({GLFW_KEY_A}, InputMode::RELEASE, WORLD_AXIS_DISPLAY);
    add_keyboard_action({GLFW_KEY_G}, InputMode::RELEASE, XY_GRID_DISPLAY);

    add_keyboard_action({GLFW_KEY_W}, InputMode::RELEASE, FACES_DISPLAY);
    add_keyboard_action({GLFW_KEY_V}, InputMode::RELEASE, VERTICES_DISPLAY);
    add_keyboard_action({GLFW_KEY_E}, InputMode::RELEASE, EDGES_DISPLAY);

    add_keyboard_action({GLFW_KEY_S}, InputMode::RELEASE, SHADING_MODE);
    add_keyboard_action({GLFW_KEY_N}, InputMode::RELEASE, INVERSE_NORMAL);
    add_keyboard_action({GLFW_KEY_M}, InputMode::RELEASE, MONO_COLOR);

    add_keyboard_action({GLFW_KEY_E, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, CYLINDER_EDGE_DISPLAY);
    add_keyboard_action({GLFW_KEY_V, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, SPHERE_VERTEX_DISPLAY);

    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD, GLFW_KEY_LEFT_CONTROL                    }, InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_6, GLFW_KEY_LEFT_CONTROL                         }, InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT, GLFW_KEY_LEFT_CONTROL               }, InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT                       }, InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD                                           }, InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_6                                                }, InputMode::HOLD, DEC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT                                      }, InputMode::HOLD, DEC_EDGES_SIZE);

    add_keyboard_action({GLFW_KEY_PAGE_UP                          }, InputMode::HOLD, INC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN                        }, InputMode::HOLD, DEC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_SHIFT  }, InputMode::HOLD, INC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_SHIFT  }, InputMode::HOLD, DEC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_ALT    }, InputMode::HOLD, INC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_ALT    }, InputMode::HOLD, DEC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, INC_LIGHT_B);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, DEC_LIGHT_B);

    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, RESET_CAMERA_AND_CP);

    /*CAMERA*/
    add_keyboard_action({GLFW_KEY_UP,   GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, FORWARD);
    add_keyboard_action({GLFW_KEY_DOWN, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, BACKWARDS);

    add_keyboard_action({GLFW_KEY_UP   }, InputMode::HOLD, UP);
    add_keyboard_action({GLFW_KEY_DOWN }, InputMode::HOLD, DOWN);
    add_keyboard_action({GLFW_KEY_LEFT }, InputMode::HOLD, LEFT);
    add_keyboard_action({GLFW_KEY_RIGHT}, InputMode::HOLD, RIGHT);

    
    add_keyboard_action({GLFW_KEY_SPACE}, InputMode::RELEASE, SWITCH_CAMERA_TYPE);
    add_keyboard_action({GLFW_KEY_O    }, InputMode::RELEASE, SWITCH_CAMERA_MODE); 

    add_keyboard_action({GLFW_KEY_A, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, SWITCH_CAMERA_CONSTRAINT_AXIS);

    add_keyboard_action({GLFW_KEY_X                     }, InputMode::HOLD, INC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R                     }, InputMode::HOLD, INC_CAMERA_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CAMERA_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT }, InputMode::HOLD, ROTATE_CAMERA);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT}, InputMode::HOLD, TRANSLATE_CAMERA);

    add_keyboard_action({GLFW_MOUSE_BUTTON_RIGHT, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, CHANGE_PIVOT_POINT);

    /*CLIPPING PLANE*/
    add_keyboard_action({GLFW_KEY_C, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, SWITCH_CLIPPING_PLANE_DISPLAY);
    add_keyboard_action({GLFW_KEY_C                   }, InputMode::RELEASE, SWITCH_CLIPPING_PLANE_MODE);

    add_keyboard_action({GLFW_KEY_A, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, SWITCH_CP_CONSTRAINT_AXIS);

    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL                     }, InputMode::HOLD, INC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL                     }, InputMode::HOLD, INC_CP_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT,   GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, ROTATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT,  GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_MIDDLE, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CP_ALONG_CAMERA_DIRECTION);

    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_TAB}, InputMode::RELEASE, RESET_CLIPPING_PLANE);

    /*ANIMATION*/
    add_keyboard_action({GLFW_KEY_F1                               }, InputMode::RELEASE, RUN_OR_STOP_ANIMATION);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_LEFT_ALT            }, InputMode::RELEASE, SAVE_KEY_FRAME);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_D, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, CLEAR_ANIMATION);

    /*===================== BIND DESCRIPTIONS ============================*/

    set_action_description({
      {"Animation", {
        {get_binding_text_from_action(SAVE_KEY_FRAME),        "Add a key frame to animation (only available in orbiter)"},
        {get_binding_text_from_action(CLEAR_ANIMATION),       "Delete animation (only available in orbiter)"},
        {get_binding_text_from_action(RUN_OR_STOP_ANIMATION), "Start/Stop animation (only available in orbiter)"},
      }},
      {"Clipping plane", {
        {get_binding_text_from_action(SWITCH_CLIPPING_PLANE_MODE),    "Switch clipping plane display mode"},
        {get_binding_text_from_action(SWITCH_CLIPPING_PLANE_DISPLAY), "Toggle clipping plane rendering on/off"},
        {get_binding_text_from_action(SWITCH_CP_CONSTRAINT_AXIS),     "Switch constraint axis for clipping plane rotation"},

        {"", ""},

        {get_binding_text_from_action(INC_CP_ROTATION_SPEED),    "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CP_ROTATION_SPEED),    "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CP_TRANSLATION_SPEED), "Increase translation speed"},
        {get_binding_text_from_action(DEC_CP_TRANSLATION_SPEED), "Decrease translation speed"},

        {"", ""},

        {get_binding_text_from_action(ROTATE_CLIPPING_PLANE),               "Rotate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CLIPPING_PLANE),            "Translate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CP_ALONG_CAMERA_DIRECTION), "Translate the clipping plane along camera direction axis when enabled"},
        {"[LCTRL+WHEEL]",                                                   "Translate the clipping plane along its normal when enabled"},

        {"", ""},

        {get_binding_text_from_action(RESET_CLIPPING_PLANE), "Reset clipping plane"},
      }},
      {"Camera", {
        {get_binding_text_from_action(FORWARD),                      "Move camera forward"},
        {get_binding_text_from_action(BACKWARDS),                    "Move camera backwards"},
        {get_binding_text_from_action(UP),                           "Move camera up"},
        {get_binding_text_from_action(DOWN),                         "Move camera down"},
        {get_binding_text_from_action(RIGHT),                        "Move camera right"},
        {get_binding_text_from_action(LEFT),                         "Move camera left"},
        {"",                                                         ""},
        {get_binding_text_from_action(SWITCH_CAMERA_MODE),           "Switch to Perspective/Orthographic view"},
        {get_binding_text_from_action(SWITCH_CAMERA_TYPE),           "Switch to Orbiter/Free-fly camera type"},
        {get_binding_text_from_action(SWITCH_CAMERA_CONSTRAINT_AXIS),"Switch constraint axis for camera rotation"},
        {"",                                                         ""},
        {get_binding_text_from_action(ROTATE_CAMERA),                "Rotate the camera"},
        {get_binding_text_from_action(TRANSLATE_CAMERA),             "Translate the camera"},
        {get_binding_text_from_action(RESET_CAMERA_AND_CP),                    "Reset camera"},
        {"",                                                         ""},
        {get_binding_text_from_action(INC_CAMERA_ROTATION_SPEED),    "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CAMERA_ROTATION_SPEED),    "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CAMERA_TRANSLATION_SPEED), "Increase translation speed"},
        {get_binding_text_from_action(DEC_CAMERA_TRANSLATION_SPEED), "Decrease translation speed"},
        {"",                                                         ""},
        {"[WHEEL]",                                                  "Zoom in/out"},
        {"[LEFT_DOUBLE_CLICK]",                                      "Aligns camera to nearest axis"},
        {"[LSHIFT+LEFT_DOUBLE_CLICK]",                               "Aligns camera to clipping plane"},
        {"[LCTRL+LEFT_DOUBLE_CLICK]",                                "Aligns clipping plane to camera"},
        {"[RIGHT_DOUBLE_CLICK]",                                     "Center the camera to the object"},
        {"[LCTRL+RIGHT_DOUBLE_CLICK]",                               "Reset camera position & orientation"},
        {get_binding_text_from_action(CHANGE_PIVOT_POINT),           "Change pivot point (center of the scene)"},
        {"[Z+WHEEL]",                                                "Increase/Decrease FOV"},
      }},
      {"Scene", {
        {get_binding_text_from_action(INC_LIGHT_ALL),         "Increase light (all colors, use shift/alt/ctrl for one rgb component)"},
        {get_binding_text_from_action(DEC_LIGHT_ALL),         "Decrease light (all colors, use shift/alt/ctrl for one rgb component)"},
        {"",                                                  ""},
        {get_binding_text_from_action(FACE_NORMALS_DISPLAY),  "Toggle face/vertex normals for normal display"},
        {get_binding_text_from_action(NORMALS_DISPLAY),       "Toggle normals display"},
        {get_binding_text_from_action(TRIANGLES_DISPLAY),     "Toggle triangles display"},
        {get_binding_text_from_action(VERTICES_DISPLAY),      "Toggle vertices display"},
        {get_binding_text_from_action(SPHERE_VERTEX_DISPLAY), "Toggle vertices display as sphere"},
        {get_binding_text_from_action(EDGES_DISPLAY),         "Toggle edges display"},
        {get_binding_text_from_action(CYLINDER_EDGE_DISPLAY), "Toggle edges display as cylinder"},
        {get_binding_text_from_action(FACES_DISPLAY),         "Toggle faces display"},
        {"",                                                  ""},
        {get_binding_text_from_action(WORLD_AXIS_DISPLAY),    "Toggle world axis display"},
        {get_binding_text_from_action(XY_GRID_DISPLAY),       "Toggle XY grid display"},
        {"",                                                  ""},
        {get_binding_text_from_action(INC_POINTS_SIZE),       "Increase size of vertices"},
        {get_binding_text_from_action(DEC_POINTS_SIZE),       "Decrease size of vertices"},
        {get_binding_text_from_action(INC_EDGES_SIZE),        "Increase size of edges"},
        {get_binding_text_from_action(DEC_EDGES_SIZE),        "Decrease size of edges"},
        {"",                                                  ""},
        {get_binding_text_from_action(MONO_COLOR),            "Toggle mono color"},
        {get_binding_text_from_action(NORMALS_MONO_COLOR),    "Toggle normals mono color"},
        {get_binding_text_from_action(INVERSE_NORMAL),        "Invert direction of normals"},
        {get_binding_text_from_action(SHADING_MODE),          "Switch between flat/Gouraud shading display"},
        {"",                                                  ""},
        {"[MIDDLE_DOUBLE_CLICK]",                             "Show entire scene"},
      }},         
      {"Window", {
        {get_binding_text_from_action(FULLSCREEN), "Switch to windowed/fullscreen mode"},
        {get_binding_text_from_action(SCREENSHOT), "Take a screenshot of the current view"},
      }},
      {"Application", {
        {get_binding_text_from_action(EXIT),                    "Exit program"},
        {"",                                                    ""},
        {get_binding_text_from_action(PRINT_APPLICATION_STATE), "Activate/Deactivate application state refreshing (FPS...)"},
      }},
    });

    print_help();
  }
} // end namespace GLFW 
} // end namespace CGAL

#endif // CGAL_BASIC_VIEWER_GLFW_IMPL_H