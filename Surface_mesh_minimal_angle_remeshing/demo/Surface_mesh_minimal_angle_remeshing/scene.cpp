#include "scene.h"
#include <fstream>
#include <QFileInfo>
#include <QOpenGLShader>
#include <QDebug>

Scene::Scene() 
  : m_frame(new ManipulatedFrame()) 
  , gl_init(false)
  , gl(NULL)
  , m_pInput(NULL)
  , m_pRemesh(NULL) {
  // 1) member data initialization
  startTimer(0);
  // 2) view option initialization
  m_view_input = false;
  m_view_remesh = true;
  m_view_facet_in_start_points = false;
  m_view_facet_in_end_points = false;
  m_view_facet_in_links = false;
  m_view_facet_out_start_points= false;
  m_view_facet_out_end_points = false;
  m_view_facet_out_links = false;
  m_view_edge_in_start_points = false;
  m_view_edge_in_end_points = false;
  m_view_edge_in_links = false;
  m_view_edge_out_start_points = false;
  m_view_edge_out_end_points = false;
  m_view_edge_out_links = false;
  m_view_vertex_in_start_points = false; 
  m_view_vertex_in_end_points = false;
  m_view_vertex_in_links = false;
  m_view_vertex_out_start_points = false;
  m_view_vertex_out_end_points = false;
  m_view_vertex_out_links = false;
  m_view_minimal_angle = true;
  m_view_polyhedron_edges = true;
  reset_draw_render_types();
  // 3) status initialization
  m_input_aabb_tree_constructed = false;
  m_links_initialized = false;
}

Scene::~Scene() {
  // member data
  delete m_frame;
  if (gl != NULL) {
    delete gl;
  }
  for (int i = 0; i < VAOs::ka_NbOfVaos; ++i) {
    vao[i].destroy();
  }
  for (int i = 0; i < VBOs::kb_NbOfVbos; ++i) {
    vbo[i].destroy();
  }
  if (m_pInput != NULL) {
    delete m_pInput;
    m_pInput = NULL;
  }
  if (m_pRemesh != NULL) {
    delete m_pRemesh;
    m_pRemesh = NULL;
  }
}

void Scene::update_bbox() {
  m_bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
    DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
  std::cout << "Computing bbox...";
  if (m_pInput == NULL) {
    std::cout << "failed (no polyhedron)." << std::endl;
    return;
  }
  if (m_pInput->empty()) {
    std::cout << "failed (empty polyhedron)." << std::endl;
    return;
  }
  Bbox bbox_input = m_pInput->get_bounding_box();
  m_bbox = m_bbox + bbox_input;
  /*if (m_pRemesh != NULL) {
  Bbox bbox_remesh = m_pRemesh->get_bounding_box();
  m_bbox = m_bbox + bbox_remesh;
  }*/
  std::cout << "Done" << std::endl;
}

void Scene::initGL() {
  gl = new QOpenGLFunctions_2_1();
  if (!gl->initializeOpenGLFunctions()) {
    qFatal("ERROR : OpenGL Functions not initialized. Check your OpenGL Verison (should be >=3.3)");
    exit(1);
  }
  compile_shaders();
  gl_init = true;
}

void Scene::compile_shaders() {
  // step 1: initialize the VAOs and VBOs
  for (int i = 0; i < VAOs::ka_NbOfVaos; ++i) {
    if (!vao[i].create()) {
      std::cerr << "VAO Creation FAILED" << std::endl;
    }
  }
  for (int i = 0; i < VBOs::kb_NbOfVbos; ++i) {
    if (!vbo[i].create()) {
      std::cerr << "VBO Creation FAILED" << std::endl;
    }
  }

  // step 2: compile the rendering_program
  const char vertex_source[] = {        // vertex source code
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 f_matrix;\n"
    "void main(void)\n"
    "{\n"
    "   gl_Position = mvp_matrix * f_matrix * vertex;\n"
    "}"
  };
  const char fragment_source[] = {      // fragment source code
    "#version 120 \n"
    "uniform highp vec4 color; \n"
    "void main(void) { \n"
    "gl_FragColor = color; \n"
    "} \n"
    "\n"
  };
  QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
  if (!vertex_shader->compileSourceCode(vertex_source)) {
    std::cerr << "Compiling vertex source FAILED" << std::endl;
  }
  QOpenGLShader *fragment_shader = new QOpenGLShader(QOpenGLShader::Fragment);
  if (!fragment_shader->compileSourceCode(fragment_source)) {
    std::cerr << "Compiling fragmentsource FAILED" << std::endl;
  }
  if (!rendering_program.addShader(vertex_shader)) {
    std::cerr << "adding vertex shader FAILED" << std::endl;
  }
  if (!rendering_program.addShader(fragment_shader)) {
    std::cerr << "adding fragment shader FAILED" << std::endl;
  }
  if (!rendering_program.link()) {
    std::cerr << "linking Program FAILED" << std::endl;
  }
  rendering_program.bind();

  // step 3: compile the rendering_program_with_light
  const char light_vertex_source[] = {        // Vertex source code
    "#version 120\n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 normals;\n"
    "attribute highp vec3 colors;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix;\n"
    "varying highp vec4 fP;\n"
    "varying highp vec3 fN;\n"
    "varying highp vec4 color;\n"
    "void main(void) {\n"
    "  color = vec4(colors, 1.0);\n"
    "  gl_Position = mvp_matrix * vertex;\n"
    "  fP = mv_matrix * vertex;\n"
    "  fN = mat3(mv_matrix)* normals;\n"
    "}\n"
  };
  const char light_fragment_source[] = {    //Fragment source code
    "#version 120\n"
    "varying highp vec4 color;\n"
    "varying highp vec4 fP;\n"
    "varying highp vec3 fN;\n"
    "uniform highp vec4 light_pos;\n"
    "uniform highp vec4 light_diff;\n"
    "uniform highp vec4 light_spec;\n"
    "uniform highp vec4 light_amb;\n"
    "uniform highp float spec_power;\n"
    "uniform int is_two_side;\n"
    "uniform bool is_selected;\n"
    "void main(void) {\n"
    "  highp vec3 L = light_pos.xyz - fP.xyz;\n"
    "  highp vec3 V = -fP.xyz;\n"
    "  highp vec3 N;\n"
    "  if (fN == highp vec3(0.0, 0.0, 0.0))\n"
    "    N = highp vec3(0.0, 0.0, 0.0);\n"
    "  else\n"
    "    N = normalize(fN);\n"
    "  L = normalize(L);\n"
    "  V = normalize(V);\n"
    "  highp vec3 R = reflect(-L, N);\n"
    "  vec4 diffuse;\n"
    "  if (is_two_side == 1)\n"
    "    diffuse = abs(dot(N, L)) * light_diff * color;\n"
    "  else\n"
    "    diffuse = max(dot(N, L), 0.0) * light_diff * color;\n"
    "  highp vec4 specular = pow(max(dot(R, V), 0.0), spec_power) * light_spec;\n"
    "  vec4 ret_color = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz, 1);\n"
    "  if (is_selected)\n"
    "    gl_FragColor = vec4(ret_color.r + 70.0 / 255.0, ret_color.g + 70.0 / 255.0, ret_color.b + 70.0 / 255.0, 1.0);\n"
    "  else\n"
    "    gl_FragColor = ret_color;\n"
    "}\n"
  };
  QOpenGLShader *light_vertex_shader = 
      new QOpenGLShader(QOpenGLShader::Vertex);
  if (!light_vertex_shader->compileSourceCode(light_vertex_source)) {
    std::cerr << "Compiling light vertex source FAILED" << std::endl;
  }
  QOpenGLShader *light_fragment_shader = 
      new QOpenGLShader(QOpenGLShader::Fragment);
  if (!light_fragment_shader->compileSourceCode(light_fragment_source)) {
    std::cerr << "Compiling light fragment source FAILED" << std::endl;
  }
  if (!rendering_program_with_light.addShader(light_vertex_shader)) {
    std::cerr << "adding light vertex shader FAILED" << std::endl;
  }
  if (!rendering_program_with_light.addShader(light_fragment_shader)) {
    std::cerr << "adding light fragment shader FAILED" << std::endl;
  }
  if (!rendering_program_with_light.link()) {
    std::cerr << "linking Program FAILED" << std::endl;
  }
  rendering_program_with_light.bind();
}

void Scene::initialize_buffers() {
  // Input facets
  vao[VAOs::ka_Input_faces].bind();
  // bind the pos_face
  vbo[VBOs::kb_Input_face_pos].bind();
  vbo[VBOs::kb_Input_face_pos].allocate(pos_input_faces.data(),
      static_cast<int>(pos_input_faces.size() * sizeof(float)));
  rendering_program_with_light.bind();
  poly_vertexLocation_with_light =
      rendering_program_with_light.attributeLocation("vertex");
  rendering_program_with_light.enableAttributeArray(
      poly_vertexLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      poly_vertexLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[VBOs::kb_Input_face_pos].release();
  // bind the pos_face_normal
  vbo[VBOs::kb_Input_face_normals].bind();
  vbo[VBOs::kb_Input_face_normals].allocate(pos_input_face_normals.data(),
      static_cast<int>(pos_input_face_normals.size() * sizeof(float)));
  rendering_program_with_light.bind();
  normalLocation_with_light =
      rendering_program_with_light.attributeLocation("normals");
  rendering_program_with_light.enableAttributeArray(
      normalLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      normalLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[VBOs::kb_Input_face_normals].release();
  // bind the pos_face_colors
  vbo[VBOs::kb_Input_face_colors].bind();
  vbo[VBOs::kb_Input_face_colors].allocate(pos_input_face_colors.data(),
    static_cast<int>(pos_input_face_colors.size() * sizeof(float)));
  rendering_program_with_light.bind();
  colorLocation_with_light =
      rendering_program_with_light.attributeLocation("colors");
  rendering_program_with_light.enableAttributeArray(colorLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      colorLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[kb_Input_face_colors].release();
  vao[VAOs::ka_Input_faces].release();

  // Remesh facets
  vao[VAOs::ka_Remesh_faces].bind();
  // bind the pos_face
  vbo[VBOs::kb_Remesh_face_pos].bind();
  vbo[VBOs::kb_Remesh_face_pos].allocate(pos_remesh_faces.data(),
      static_cast<int>(pos_remesh_faces.size() * sizeof(float)));
  rendering_program_with_light.bind();
  poly_vertexLocation_with_light =
      rendering_program_with_light.attributeLocation("vertex");
  rendering_program_with_light.enableAttributeArray(
      poly_vertexLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      poly_vertexLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[VBOs::kb_Remesh_face_pos].release();
  // bind the pos_face_normal
  vbo[VBOs::kb_Remesh_face_normals].bind();
  vbo[VBOs::kb_Remesh_face_normals].allocate(pos_remesh_face_normals.data(),
      static_cast<int>(pos_remesh_face_normals.size() * sizeof(float)));
  rendering_program_with_light.bind();
  normalLocation_with_light =
      rendering_program_with_light.attributeLocation("normals");
  rendering_program_with_light.enableAttributeArray(normalLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      normalLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[VBOs::kb_Remesh_face_normals].release();
  // bind the pos_ifi_colors
  vbo[VBOs::kb_Remesh_face_colors].bind();
  vbo[VBOs::kb_Remesh_face_colors].allocate(pos_remesh_face_colors.data(),
      static_cast<int>(pos_remesh_face_colors.size() * sizeof(float)));
  rendering_program_with_light.bind();
  colorLocation_with_light =
      rendering_program_with_light.attributeLocation("colors");
  rendering_program_with_light.enableAttributeArray(colorLocation_with_light);
  rendering_program_with_light.setAttributeBuffer(
      colorLocation_with_light, GL_FLOAT, 0, 3);
  rendering_program_with_light.release();
  vbo[kb_Remesh_face_colors].release();
  vao[VAOs::ka_Remesh_faces].release();

  // Input edges
  vao[VAOs::ka_Input_boundaries].bind();
  vbo[VBOs::kb_Input_boundary_pos].bind();
  vbo[VBOs::kb_Input_boundary_pos].allocate(pos_input_boundaries.data(),
      static_cast<int>(pos_input_boundaries.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Input_boundary_pos].release();
  vao[VAOs::ka_Input_boundaries].release();

  // Remesh edges
  vao[VAOs::ka_Remesh_boundaries].bind();
  vbo[VBOs::kb_Remesh_boundary_pos].bind();
  vbo[VBOs::kb_Remesh_boundary_pos].allocate(pos_remesh_boundaries.data(),
      static_cast<int>(pos_remesh_boundaries.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Remesh_boundary_pos].release();
  vao[VAOs::ka_Remesh_boundaries].release();

  // Input samples
  vao[VAOs::ka_Input_samples].bind();
  vbo[VBOs::kb_Input_sample_pos].bind();
  vbo[VBOs::kb_Input_sample_pos].allocate(pos_input_samples.data(),
      static_cast<int>(pos_input_samples.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Input_sample_pos].release();
  rendering_program.release();
  vao[VAOs::ka_Input_samples].release();

  // Remesh samples
  vao[VAOs::ka_Remesh_samples].bind();
  vbo[VBOs::kb_Remesh_sample_pos].bind();
  vbo[VBOs::kb_Remesh_sample_pos].allocate(pos_remesh_samples.data(),
      static_cast<int>(pos_remesh_samples.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Remesh_sample_pos].release();
  rendering_program.release();
  vao[VAOs::ka_Remesh_samples].release();

  // Input normal edges
  vao[VAOs::ka_Input_normal_edges].bind();
  vbo[VBOs::kb_Input_normal_edge_pos].bind();
  vbo[VBOs::kb_Input_normal_edge_pos].allocate(pos_input_normal_edges.data(),
      static_cast<int>(pos_input_normal_edges.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Input_normal_edge_pos].release();
  vao[VAOs::ka_Input_normal_edges].release();

  // Remesh normal edges
  vao[VAOs::ka_Remesh_normal_edges].bind();
  vbo[VBOs::kb_Remesh_normal_edge_pos].bind();
  vbo[VBOs::kb_Remesh_normal_edge_pos].allocate(pos_remesh_normal_edges.data(),
      static_cast<int>(pos_remesh_normal_edges.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Remesh_normal_edge_pos].release();
  vao[VAOs::ka_Remesh_normal_edges].release();

  // Input min radian edges
  vao[VAOs::ka_Input_special_edges].bind();
  vbo[VBOs::kb_Input_special_edge_pos].bind();
  vbo[VBOs::kb_Input_special_edge_pos].allocate(pos_input_special_edges.data(),
      static_cast<int>(pos_input_special_edges.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Input_special_edge_pos].release();
  vao[VAOs::ka_Input_special_edges].release();

  // Remesh min radian edges
  vao[VAOs::ka_Remesh_special_edges].bind();
  vbo[VBOs::kb_Remesh_special_edge_pos].bind();
  vbo[VBOs::kb_Remesh_special_edge_pos].allocate(
      pos_remesh_special_edges.data(),
      static_cast<int>(pos_remesh_special_edges.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  rendering_program.release();
  vbo[VBOs::kb_Remesh_special_edge_pos].release();
  vao[VAOs::ka_Remesh_special_edges].release();

  // Facet in start points
  vao[VAOs::ka_Facet_in_start].bind();
  vbo[VBOs::kb_Facet_in_start_points].bind();
  vbo[VBOs::kb_Facet_in_start_points].allocate(
      pos_facet_in_start_points.data(),
      static_cast<int>(pos_facet_in_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Facet_in_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_in_start].release();

  // Facet in end points
  vao[VAOs::ka_Facet_in_end].bind();
  vbo[VBOs::kb_Facet_in_end_points].bind();
  vbo[VBOs::kb_Facet_in_end_points].allocate(pos_facet_in_end_points.data(),
      static_cast<int>(pos_facet_in_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Facet_in_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_in_end].release();

  // Facet in links
  vao[VAOs::ka_Facet_in_links].bind();
  vbo[VBOs::kb_Facet_in_link_lines].bind();
  vbo[VBOs::kb_Facet_in_link_lines].allocate(pos_facet_in_links.data(),
      static_cast<int>(pos_facet_in_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Facet_in_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_in_links].release();

  // Facet out start points
  vao[VAOs::ka_Facet_out_start].bind();
  vbo[VBOs::kb_Facet_out_start_points].bind();
  vbo[VBOs::kb_Facet_out_start_points].allocate(
      pos_facet_out_start_points.data(),
      static_cast<int>(pos_facet_out_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Facet_out_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_out_start].release();

  // Facet out end points
  vao[VAOs::ka_Facet_out_end].bind();
  vbo[VBOs::kb_Facet_out_end_points].bind();
  vbo[VBOs::kb_Facet_out_end_points].allocate(pos_facet_out_end_points.data(),
      static_cast<int>(pos_facet_out_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Facet_out_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_out_end].release();

  // Facet out links
  vao[VAOs::ka_Facet_out_links].bind();
  vbo[VBOs::kb_Facet_out_link_lines].bind();
  vbo[VBOs::kb_Facet_out_link_lines].allocate(pos_facet_out_links.data(),
      static_cast<int>(pos_facet_out_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Facet_out_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Facet_out_links].release();

  // Edge in start points
  vao[VAOs::ka_Edge_in_start].bind();
  vbo[VBOs::kb_Edge_in_start_points].bind();
  vbo[VBOs::kb_Edge_in_start_points].allocate(pos_edge_in_start_points.data(),
      static_cast<int>(pos_edge_in_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Edge_in_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_in_start].release();

  // Edge in end points
  vao[VAOs::ka_Edge_in_end].bind();
  vbo[VBOs::kb_Edge_in_end_points].bind();
  vbo[VBOs::kb_Edge_in_end_points].allocate(pos_edge_in_end_points.data(),
      static_cast<int>(pos_edge_in_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Edge_in_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_in_end].release();

  // Edge in links
  vao[VAOs::ka_Edge_in_links].bind();
  vbo[VBOs::kb_Edge_in_link_lines].bind();
  vbo[VBOs::kb_Edge_in_link_lines].allocate(pos_edge_in_links.data(),
      static_cast<int>(pos_edge_in_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Edge_in_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_in_links].release();

  // Edge out start points
  vao[VAOs::ka_Edge_out_start].bind();
  vbo[VBOs::kb_Edge_out_start_points].bind();
  vbo[VBOs::kb_Edge_out_start_points].allocate(
      pos_edge_out_start_points.data(),
      static_cast<int>(pos_edge_out_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Edge_out_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_out_start].release();

  // Edge out end points
  vao[VAOs::ka_Edge_out_end].bind();
  vbo[VBOs::kb_Edge_out_end_points].bind();
  vbo[VBOs::kb_Edge_out_end_points].allocate(pos_edge_out_end_points.data(),
      static_cast<int>(pos_edge_out_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Edge_out_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_out_end].release();

  // Edge out links
  vao[VAOs::ka_Edge_out_links].bind();
  vbo[VBOs::kb_Edge_out_link_lines].bind();
  vbo[VBOs::kb_Edge_out_link_lines].allocate(pos_edge_out_links.data(),
      static_cast<int>(pos_edge_out_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Edge_out_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Edge_out_links].release();

  // Vertex in start points
  vao[VAOs::ka_Vertex_in_start].bind();
  vbo[VBOs::kb_Vertex_in_start_points].bind();
  vbo[VBOs::kb_Vertex_in_start_points].allocate(
      pos_vertex_in_start_points.data(),
      static_cast<int>(pos_vertex_in_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Vertex_in_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_in_start].release();

  // Vertex in end points
  vao[VAOs::ka_Vertex_in_end].bind();
  vbo[VBOs::kb_Vertex_in_end_points].bind();
  vbo[VBOs::kb_Vertex_in_end_points].allocate(pos_vertex_in_end_points.data(),
      static_cast<int>(pos_vertex_in_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Vertex_in_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_in_end].release();

  // Vertex in links
  vao[VAOs::ka_Vertex_in_links].bind();
  vbo[VBOs::kb_Vertex_in_link_lines].bind();
  vbo[VBOs::kb_Vertex_in_link_lines].allocate(pos_vertex_in_links.data(),
      static_cast<int>(pos_vertex_in_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Vertex_in_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_in_links].release();

  // Vertex out start points
  vao[VAOs::ka_Vertex_out_start].bind();
  vbo[VBOs::kb_Vertex_out_start_points].bind();
  vbo[VBOs::kb_Vertex_out_start_points].allocate(
      pos_vertex_out_start_points.data(),
      static_cast<int>(pos_vertex_out_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Vertex_out_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_out_start].release();

  // Vertex out end points
  vao[VAOs::ka_Vertex_out_end].bind();
  vbo[VBOs::kb_Vertex_out_end_points].bind();
  vbo[VBOs::kb_Vertex_out_end_points].allocate( 
      pos_vertex_out_end_points.data(),
      static_cast<int>(pos_vertex_out_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Vertex_out_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_out_end].release();

  // Vertex out links
  vao[VAOs::ka_Vertex_out_links].bind();
  vbo[VBOs::kb_Vertex_out_link_lines].bind();
  vbo[VBOs::kb_Vertex_out_link_lines].allocate(pos_vertex_out_links.data(),
      static_cast<int>(pos_vertex_out_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Vertex_out_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Vertex_out_links].release();

  are_buffers_initialized = true;
}

void Scene::attrib_buffers(CGAL::QGLViewer *viewer) {
  QMatrix4x4 mvMatrix, mvpMatrix;
  double mvMat[16], mvpMat[16];
  viewer->camera()->getModelViewMatrix(mvMat);
  viewer->camera()->getModelViewProjectionMatrix(mvpMat);
  for (int i = 0; i < 16; ++i) {
    mvMatrix.data()[i] = (float)mvMat[i];
    mvpMatrix.data()[i] = (float)mvpMat[i];
  }

  rendering_program.bind();
  mvpLocation = rendering_program.uniformLocation("mvp_matrix");
  fLocation = rendering_program.uniformLocation("f_matrix");
  colorLocation = rendering_program.uniformLocation("color");
  rendering_program.setUniformValue(mvpLocation, mvpMatrix);
  rendering_program.release();

  QVector4D position(0.0f, 0.0f, 1.0f, 1.0f);
  //QVector4D diffuse(0.8f, 0.8f, 0.8f, 1.0f);    // backup
  //QVector4D specular(0.2f, 0.2f, 0.2f, 1.0f);
  //QVector4D ambient(0.8f, 0.8f, 0.8f, 1.0f);
  //float spec_power = 51.8f;
  QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
  QVector4D specular(0.05f, 0.05f, 0.05f, 1.0f);
  QVector4D ambient(0.4f, 0.4f, 0.4f, 1.0f);
  float spec_power = 20.0f;
  GLint is_both_sides = 0;
  rendering_program_with_light.bind();
  rendering_program_with_light.setUniformValue("light_pos", position);
  rendering_program_with_light.setUniformValue("light_diff", diffuse);
  rendering_program_with_light.setUniformValue("light_spec", specular);
  rendering_program_with_light.setUniformValue("light_amb", ambient);
  rendering_program_with_light.setUniformValue("spec_power", spec_power);
  rendering_program_with_light.setUniformValue("is_two_side", is_both_sides);
  rendering_program_with_light.setUniformValue("mv_matrix", mvMatrix);
  rendering_program_with_light.setUniformValue("mvp_matrix", mvpMatrix);
  rendering_program_with_light.release();
}

void Scene::draw(CGAL::QGLViewer *viewer) {
  if (!gl_init) {
    initGL();
  }
  if (!are_buffers_initialized) {
    initialize_buffers();
  }
  gl->glEnable(GL_DEPTH_TEST);
  QColor color;
  QMatrix4x4 fMatrix;
  fMatrix.setToIdentity();
  // Input
  if (m_view_input) {
    if (pos_input_faces.size() > 0) {
      gl->glEnable(GL_LIGHTING);
      gl->glEnable(GL_POLYGON_OFFSET_FILL);
      gl->glPolygonOffset(1.0f, 1.0f);
      vao[VAOs::ka_Input_faces].bind();
      attrib_buffers(viewer);
      rendering_program_with_light.bind();
      rendering_program_with_light.setUniformValue("is_selected", false);
      rendering_program_with_light.setUniformValue("f_matrix", fMatrix);
      gl->glDrawArrays(GL_TRIANGLES, 0,
        static_cast<GLsizei>(pos_input_faces.size() / 3));
      rendering_program_with_light.release();
      vao[VAOs::ka_Input_faces].release();
      gl->glDisable(GL_POLYGON_OFFSET_FILL);
    }
    if (pos_input_boundaries.size() > 0) {
      gl->glDisable(GL_LIGHTING);
      gl->glLineWidth(0.5f);
      vao[VAOs::ka_Input_boundaries].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 200, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_input_boundaries.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Input_boundaries].release();
    }
    if (pos_input_samples.size() > 0) {
      gl->glDisable(GL_LIGHTING);
      gl->glPointSize(3.0f);
      vao[VAOs::ka_Input_samples].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 255, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_input_samples.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Input_samples].release();
    }
    if (pos_input_normal_edges.size() > 0 &&
      (m_render_type == k_classifications || m_view_polyhedron_edges)) {
      gl->glDisable(GL_LIGHTING);
      if (m_render_type == k_classifications || pos_input_boundaries.empty()) {
        gl->glLineWidth(1.0f);
      }
      else {
        gl->glLineWidth(2.0f);
      }
      vao[VAOs::ka_Input_normal_edges].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      m_render_type == k_classifications ? color.setRgb(0, 0, 255) :
        color.setRgb(0, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_input_normal_edges.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Input_normal_edges].release();
    }
    if (pos_input_special_edges.size() > 0 &&
      (m_render_type == k_classifications ||
      (m_view_minimal_angle && m_render_type == k_plain_facets))) {
      gl->glDisable(GL_LIGHTING);
      gl->glLineWidth(2.0f);
      vao[VAOs::ka_Input_special_edges].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(255, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_input_special_edges.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Input_special_edges].release();
    }
  }
  // Remesh
  if (m_view_remesh) {
    if (pos_remesh_faces.size() > 0) {
      gl->glEnable(GL_LIGHTING);
      gl->glEnable(GL_POLYGON_OFFSET_FILL);
      gl->glPolygonOffset(1.0f, 1.0f);
      vao[VAOs::ka_Remesh_faces].bind();
      attrib_buffers(viewer);
      rendering_program_with_light.bind();
      rendering_program_with_light.setUniformValue("is_selected", false);
      rendering_program_with_light.setUniformValue("f_matrix", fMatrix);
      gl->glDrawArrays(GL_TRIANGLES, 0,
        static_cast<GLsizei>(pos_remesh_faces.size() / 3));
      rendering_program_with_light.release();
      vao[VAOs::ka_Remesh_faces].release();
      gl->glDisable(GL_POLYGON_OFFSET_FILL);
    }
    if (pos_remesh_boundaries.size() > 0) {
      gl->glDisable(GL_LIGHTING);
      gl->glLineWidth(0.5f);
      vao[VAOs::ka_Remesh_boundaries].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 200, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_remesh_boundaries.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Remesh_boundaries].release();
    }
    if (pos_remesh_samples.size() > 0) {
      gl->glDisable(GL_LIGHTING);
      gl->glPointSize(3.0f);
      vao[VAOs::ka_Remesh_samples].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 255, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_remesh_samples.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Remesh_samples].release();
    }
    if (pos_remesh_normal_edges.size() > 0 &&
      (m_render_type == k_classifications || m_view_polyhedron_edges)) {
      gl->glDisable(GL_LIGHTING);
      if (m_render_type == k_classifications || 
          pos_remesh_boundaries.empty()) {
        gl->glLineWidth(1.0f);
      }
      else {
        gl->glLineWidth(2.0f);
      }
      vao[VAOs::ka_Remesh_normal_edges].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      m_render_type == k_classifications ? color.setRgb(0, 0, 255) :
        color.setRgb(0, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_remesh_normal_edges.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Remesh_normal_edges].release();
    }
    if (pos_remesh_special_edges.size() > 0 &&
      (m_render_type == k_classifications ||
      (m_view_minimal_angle && m_render_type == k_plain_facets))) {
      gl->glDisable(GL_LIGHTING);
      gl->glLineWidth(2.0f);
      vao[VAOs::ka_Remesh_special_edges].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(255, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_remesh_special_edges.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Remesh_special_edges].release();
    }
  }
  // Samples and links
  if (m_draw_type == DrawType::k_polyhedron) {
    gl->glDisable(GL_LIGHTING);
    gl->glPointSize(3.0f);
    gl->glLineWidth(1.0f);
    // facet in start points
    if (m_view_facet_in_start_points && pos_facet_in_start_points.size() > 0) {
      vao[VAOs::ka_Facet_in_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 250, 250);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_facet_in_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_in_start].release();
    }
    // facet in links
    if (m_view_facet_in_links && pos_facet_in_links.size() > 0) {
      vao[VAOs::ka_Facet_in_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 200, 200);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_facet_in_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_in_links].release();
    }
    // facet in end points
    if (m_view_facet_in_end_points && pos_facet_in_end_points.size() > 0) {
      vao[VAOs::ka_Facet_in_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 150, 150);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_facet_in_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_in_end].release();
    }
    // facet out start points
    if (m_view_facet_out_start_points && 
        pos_facet_out_start_points.size() > 0) {
      vao[VAOs::ka_Facet_out_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(250, 0, 250);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_facet_out_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_out_start].release();
    }
    // facet out links
    if (m_view_facet_out_links && pos_facet_out_links.size() > 0) {
      vao[VAOs::ka_Facet_out_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(200, 0, 200);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_facet_out_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_out_links].release();
    }
    // facet out end points
    if (m_view_facet_out_end_points && pos_facet_out_end_points.size() > 0) {
      vao[VAOs::ka_Facet_out_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(150, 0, 150);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_facet_out_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Facet_out_end].release();
    }
    // edge in start points
    if (m_view_edge_in_start_points && pos_edge_in_start_points.size() > 0) {
      vao[VAOs::ka_Edge_in_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 0, 250);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_edge_in_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_in_start].release();
    }
    // edge in links
    if (m_view_edge_in_links && pos_edge_in_links.size() > 0) {
      vao[VAOs::ka_Edge_in_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 0, 200);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_edge_in_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_in_links].release();
    }
    // edge in end points
    if (m_view_edge_in_end_points && pos_edge_in_end_points.size() > 0) {
      vao[VAOs::ka_Edge_in_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 0, 150);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_edge_in_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_in_end].release();
    }
    // edge out start points
    if (m_view_edge_out_start_points && pos_edge_out_start_points.size() > 0) {
      vao[VAOs::ka_Edge_out_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(250, 250, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_edge_out_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_out_start].release();
    }
    // edge out links
    if (m_view_edge_out_links && pos_edge_out_links.size() > 0) {
      vao[VAOs::ka_Edge_out_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(200, 200, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_edge_out_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_out_links].release();
    }
    // edge out end points
    if (m_view_edge_out_end_points && pos_edge_out_end_points.size() > 0) {
      vao[VAOs::ka_Edge_out_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(150, 150, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_edge_out_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Edge_out_end].release();
    }
    // vertex in start points
    if (m_view_vertex_in_start_points &&
        pos_vertex_in_start_points.size() > 0) {
      vao[VAOs::ka_Vertex_in_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 250, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_vertex_in_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_in_start].release();
    }
    // vertex in links
    if (m_view_vertex_in_links && pos_vertex_in_links.size() > 0) {
      vao[VAOs::ka_Vertex_in_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 200, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_vertex_in_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_in_links].release();
    }
    // vertex in end points
    if (m_view_vertex_in_end_points && pos_vertex_in_end_points.size() > 0) {
      vao[VAOs::ka_Vertex_in_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 150, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_vertex_in_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_in_end].release();
    }
    // vertex out start points
    if (m_view_vertex_out_start_points &&
        pos_vertex_out_start_points.size() > 0) {
      vao[VAOs::ka_Vertex_out_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(250, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_vertex_out_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_out_start].release();
    }
    // vertex out links
    if (m_view_vertex_out_links && pos_vertex_out_links.size() > 0) {
      vao[VAOs::ka_Vertex_out_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(200, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_vertex_out_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_out_links].release();
    }
    // vertex out end points
    if (m_view_vertex_out_end_points && pos_vertex_out_end_points.size() > 0) {
      vao[VAOs::ka_Vertex_out_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(150, 0, 0);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_vertex_out_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Vertex_out_end].release();
    }
  }
}

void Scene::changed() {
  compute_elements();
  are_buffers_initialized = false;
}

void Scene::compute_elements() {
  if (m_view_input && m_pInput != NULL) {
    // step 1: compute the facets
    remesh.compute_facets(m_draw_type, m_render_type, m_bbox, true, m_pInput,
        &pos_input_faces, &pos_input_face_normals, &pos_input_face_colors,
        &pos_input_boundaries, &pos_input_samples);
    // step 2: render the edges
    if (m_render_type == k_classifications) {
      m_pInput->compute_classified_edges(&pos_input_normal_edges, 
                                         &pos_input_special_edges);
    }
    else {
      if (m_view_polyhedron_edges) {
        m_pInput->compute_edges(&pos_input_normal_edges);
      }
      if (m_view_minimal_angle && m_render_type == k_plain_facets) {
        m_pInput->compute_min_radian_edges(&pos_input_special_edges);
      }
    }
  }
  if (m_view_remesh && m_pRemesh != NULL) {
    // step 1: render the facets
    remesh.compute_facets(m_draw_type, m_render_type, m_bbox, false, m_pRemesh,
        &pos_remesh_faces, &pos_remesh_face_normals, &pos_remesh_face_colors,
        &pos_remesh_boundaries, &pos_remesh_samples);
    // step 2: render the edges
    if (m_render_type == k_classifications) {
      m_pRemesh->compute_classified_edges(&pos_remesh_normal_edges,
                                          &pos_remesh_special_edges);
    }
    else {
      if (m_view_polyhedron_edges) {
        m_pRemesh->compute_edges(&pos_remesh_normal_edges);
      }
      if (m_view_minimal_angle && m_render_type == k_plain_facets) {
        m_pRemesh->compute_min_radian_edges(&pos_remesh_special_edges);
      }
    }
  }
  // step 3: compute samples and links
  if (m_pInput != NULL && m_pRemesh != NULL && m_draw_type == k_polyhedron) {
    // facet in links
    if (m_view_facet_in_start_points) {
      m_pInput->compute_facet_start_points(&pos_facet_in_start_points);
    }
    if (m_view_facet_in_end_points) {
      m_pInput->compute_facet_end_points(&pos_facet_in_end_points);
    }
    if (m_view_facet_in_links) {
      m_pInput->compute_facet_links(&pos_facet_in_links);
    }
    // edge in links
    if (m_view_edge_in_start_points) {
      m_pInput->compute_edge_start_points(&pos_edge_in_start_points);
    }
    if (m_view_edge_in_end_points) {
      m_pInput->compute_edge_end_points(&pos_edge_in_end_points);
    }
    if (m_view_edge_in_links) {
      m_pInput->compute_edge_links(&pos_edge_in_links);
    }
    // vertex in links
    if (m_view_vertex_in_start_points) {
      m_pInput->compute_vertex_start_points(&pos_vertex_in_start_points);
    }
    if (m_view_vertex_in_end_points) {
      m_pInput->compute_vertex_end_points(&pos_vertex_in_end_points);
    }
    if (m_view_vertex_in_links) {
      m_pInput->compute_vertex_links(&pos_vertex_in_links);
    }
    // facet out links
    if (m_view_facet_out_start_points) {
      m_pRemesh->compute_facet_start_points(&pos_facet_out_start_points);
    }
    if (m_view_facet_out_end_points) {
      m_pRemesh->compute_facet_end_points(&pos_facet_out_end_points);
    }
    if (m_view_facet_out_links) {
      m_pRemesh->compute_facet_links(&pos_facet_out_links);
    }
    // edge out links
    if (m_view_edge_out_start_points) {
      m_pRemesh->compute_edge_start_points(&pos_edge_out_start_points);
    }
    if (m_view_edge_out_end_points) {
      m_pRemesh->compute_edge_end_points(&pos_edge_out_end_points);
    }
    if (m_view_edge_out_links) {
      m_pRemesh->compute_edge_links(&pos_edge_out_links);
    }
    // vertex out links
    if (m_view_vertex_out_start_points) {
      m_pRemesh->compute_vertex_start_points(&pos_vertex_out_start_points);
    }
    if (m_view_vertex_out_end_points) {
      m_pRemesh->compute_vertex_end_points(&pos_vertex_out_end_points);
    }
    if (m_view_vertex_out_links) {
      m_pRemesh->compute_vertex_links(&pos_vertex_out_links);
    }
  }
}

bool Scene::open(QString file_name) {
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"\n").arg(file_name);
  cerr.flush();

  QFileInfo file_info(file_name);
  std::ifstream in(file_name.toUtf8());
  if (!file_info.isFile() || !file_info.isReadable() || !in) {
    std::cerr << "unable to open file" << std::endl;
    return false;
  }

  if (m_pInput != NULL) {
    delete m_pInput;
  }
  if (m_pRemesh != NULL) {
    delete m_pRemesh;
  }

  m_pInput = new Polyhedron;
  in >> *m_pInput;
  if (!in) {
    std::cerr << "invalid OFF file" << std::endl;
    delete m_pInput;
    m_pInput = NULL;
    return false;
  }
  in.close();

  m_pInput->normalize(1.0);
  m_pInput->calculate_normals("Input");
  remesh.calculate_feature_intensities("Input", m_pInput);
  update_bbox();
  m_pRemesh = new Polyhedron(*m_pInput);
  m_input_aabb_tree_constructed = false;
  m_links_initialized = false;
  remesh.initialize_private_data();

  reset_draw_render_types();
  changed();

  return true;
}

bool Scene::open_input(QString file_name) {
  QTextStream cerr(stderr);
  cerr << QString("Opening input file \"%1\"\n").arg(file_name);
  cerr.flush();

  QFileInfo fileinfo(file_name);
  std::ifstream in(file_name.toUtf8());
  if (!in || !fileinfo.isFile() || !fileinfo.isReadable()) {
    std::cerr << "unable to open input file" << std::endl;
    return false;
  }

  if (m_pInput != NULL) {
    delete m_pInput;
  }
  m_pInput = new Polyhedron;
  in >> *m_pInput;
  if (!in) {
    std::cerr << "invalid OFF file" << std::endl;
    delete m_pInput;
    m_pInput = NULL;
    return false;
  }
  in.close();

  m_pInput->normalize(1.0);
  m_pInput->calculate_normals("Input");
  remesh.calculate_feature_intensities("Input", m_pInput);
  update_bbox();
  m_input_aabb_tree_constructed = false;

  if (m_pRemesh == NULL) {
    m_pRemesh = new Polyhedron(*m_pInput);
  }
  remesh.clear_links(m_pInput, m_pRemesh);
  m_links_initialized = false;
  remesh.initialize_private_data();

  reset_draw_render_types();
  m_view_input = true;
  m_view_remesh = false;
  changed();

  return true;
}

bool Scene::open_remesh(QString file_name) {
  QTextStream cerr(stderr);
  cerr << QString("Opening remesh file \"%1\"\n").arg(file_name);
  cerr.flush();

  QFileInfo fileinfo(file_name);
  std::ifstream in(file_name.toUtf8());
  if (!in || !fileinfo.isFile() || !fileinfo.isReadable()) {
    std::cerr << "unable to open remesh file" << std::endl;
    return false;
  }

  if (m_pRemesh != NULL) {
    delete m_pRemesh;
  }
  m_pRemesh = new Polyhedron;
  in >> *m_pRemesh;
  if (!in) {
    std::cerr << "invalid OFF file" << std::endl;
    delete m_pRemesh;
    m_pRemesh = NULL;
    return false;
  }
  in.close();

  m_pRemesh->calculate_normals("Remesh");
  remesh.calculate_feature_intensities("Remesh", m_pRemesh);
  if (m_pInput == NULL) {
    m_pRemesh->normalize(1.0);
    m_pInput = new Polyhedron(*m_pRemesh);
    m_input_aabb_tree_constructed = false;
    update_bbox();
  }
  remesh.clear_links(m_pInput, m_pRemesh);
  m_links_initialized = false;
  remesh.initialize_private_data();

  reset_draw_render_types();
  changed();

  return true;
}

void Scene::save_remesh_as(QString file_name) {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_pRemesh->save_as(file_name.toStdString());
  }
}

int Scene::get_optimize_type_index(OptimizeType ot) const {
  switch (ot) {
  case OptimizeType::k_none:
    return 0;
  case OptimizeType::k_input_to_remesh:
    return 1;
  case OptimizeType::k_remesh_to_input:
    return 2;
  case OptimizeType::k_both:
    return 3;
  default:
    return -1;    // invalid case
  }
}

OptimizeType Scene::get_optimize_type(int index) const {
  switch (index) {
  case 0:
    return OptimizeType::k_none;
  case 1:
    return OptimizeType::k_input_to_remesh;
  case 2:
    return OptimizeType::k_remesh_to_input;
  case 3:
    return OptimizeType::k_both;
  default:
    return OptimizeType::k_none;
  }
}

void Scene::toggle_view_input() { 
  m_view_input = !m_view_input;
  if (m_view_input) {
    changed();
  }
}

void Scene::toggle_view_remesh() { 
  m_view_remesh = !m_view_remesh; 
  if (m_view_remesh) {
    changed();
  }
}

void Scene::toggle_view_input_remesh() {
  m_view_remesh = !m_view_remesh;
  m_view_input = !m_view_remesh;
  if (m_view_input || m_view_remesh) {
    changed();
  }
}

void Scene::toggle_view_minimal_angle() {
  m_view_minimal_angle = !m_view_minimal_angle;
  if (m_view_minimal_angle) {
    changed();
  }
}

void Scene::toggle_view_polyhedron_edges() {
  m_view_polyhedron_edges = !m_view_polyhedron_edges;
  if (m_view_polyhedron_edges) {
    changed();
  }
}

void Scene::toggle_view_polyhedron_facets() {
  reset_draw_render_types();
  changed();
}

void Scene::toggle_view_interpolated_feature_intensities() {
  set_draw_render_types(k_polyhedron, k_ifi_facets);
  changed();
}

void Scene::toggle_view_facet_errors() {
  set_draw_render_types(k_polyhedron, k_mr_facets);
  changed();
}

void Scene::toggle_view_element_classifications() {
  set_draw_render_types(k_vertex_voronoi, k_classifications);
  changed();
}

void Scene::toggle_view_gaussian_curvatures() {
  set_draw_render_types(k_vertex_voronoi, k_gaussian_curvature);
  changed();
}

void Scene::toggle_view_maximal_normal_dihedrals() {
  set_draw_render_types(k_vertex_voronoi, k_maximal_halfedge_dihedral);
  changed();
}

void Scene::toggle_view_normal_dihedrals() {
  set_draw_render_types(k_edge_voronoi, k_normal_dihedral);
  changed();
}

void Scene::toggle_view_facet_in_start_points() {
  m_view_facet_in_start_points = !m_view_facet_in_start_points;
  if (m_view_facet_in_start_points) {
    changed();
  }
}

void Scene::toggle_view_facet_in_end_points() {
  m_view_facet_in_end_points = !m_view_facet_in_end_points;
  if (m_view_facet_in_end_points) {
    changed();
  }
}

void Scene::toggle_view_facet_in_links() {
  m_view_facet_in_links = !m_view_facet_in_links;
  if (m_view_facet_in_links) {
    changed();
  }
  changed();
}

void Scene::toggle_view_facet_out_start_points() {
  m_view_facet_out_start_points = !m_view_facet_out_start_points;
  if (m_view_facet_out_start_points) {
    changed();
  }
}

void Scene::toggle_view_facet_out_end_points() {
  m_view_facet_out_end_points = !m_view_facet_out_end_points;
  if (m_view_facet_out_end_points) {
    changed();
  }
}

void Scene::toggle_view_facet_out_links() {
  m_view_facet_out_links = !m_view_facet_out_links;
  if (m_view_facet_out_links) {
    changed();
  }
}

void Scene::toggle_view_edge_in_start_points() {
  m_view_edge_in_start_points = !m_view_edge_in_start_points;
  if (m_view_edge_in_start_points) {
    changed();
  }
}

void Scene::toggle_view_edge_in_end_points() {
  m_view_edge_in_end_points = !m_view_edge_in_end_points;
  if (m_view_edge_in_end_points) {
    changed();
  }
}

void Scene::toggle_view_edge_in_links() {
  m_view_edge_in_links = !m_view_edge_in_links;
  if (m_view_edge_in_links) {
    changed();
  }
}

void Scene::toggle_view_edge_out_start_points() {
  m_view_edge_out_start_points = !m_view_edge_out_start_points;
  if (m_view_edge_out_start_points) {
    changed();
  }
}

void Scene::toggle_view_edge_out_end_points() {
  m_view_edge_out_end_points = !m_view_edge_out_end_points;
  if (m_view_edge_out_end_points) {
    changed();
  }
}

void Scene::toggle_view_edge_out_links() {
  m_view_edge_out_links = !m_view_edge_out_links;
  if (m_view_edge_out_links) {
    changed();
  }
}

void Scene::toggle_view_vertex_in_start_points() {
  m_view_vertex_in_start_points = !m_view_vertex_in_start_points;
  if (m_view_vertex_in_start_points) {
    changed();
  }
}

void Scene::toggle_view_vertex_in_end_points() {
  m_view_vertex_in_end_points = !m_view_vertex_in_end_points;
  if (m_view_vertex_in_end_points) {
    changed();
  }
}

void Scene::toggle_view_vertex_in_links() {
  m_view_vertex_in_links = !m_view_vertex_in_links;
  if (m_view_vertex_in_links) {
    changed();
  }
}

void Scene::toggle_view_vertex_out_start_points() {
  m_view_vertex_out_start_points = !m_view_vertex_out_start_points;
  if (m_view_vertex_out_start_points) {
    changed();
  }
}

void Scene::toggle_view_vertex_out_end_points() {
  m_view_vertex_out_end_points = !m_view_vertex_out_end_points;
  if (m_view_vertex_out_end_points) {
    changed();
  }
}

void Scene::toggle_view_vertex_out_links() {
  m_view_vertex_out_links = !m_view_vertex_out_links;
  if (m_view_vertex_out_links) {
    changed();
  }
}

void Scene::toggle_view_all_sample_feature_intensities() {  // all samples
  set_draw_render_types(k_all_voronoi, k_feature_intensity);
  changed();
}

void Scene::toggle_view_all_sample_capacities() {
  set_draw_render_types(k_all_voronoi, k_capacity);
  changed();
}

void Scene::toggle_view_all_sample_weights() {
  set_draw_render_types(k_all_voronoi, k_weight);
  changed();
}

void Scene::toggle_view_vertex_feature_intensities() {      // vertex samples
  set_draw_render_types(k_vertex_voronoi, k_feature_intensity);
  changed();
}

void Scene::toggle_view_vertex_capacities() {
  set_draw_render_types(k_vertex_voronoi, k_capacity);
  changed();
}

void Scene::toggle_view_vertex_weights() {
  set_draw_render_types(k_vertex_voronoi, k_weight);
  changed();
}

void Scene::toggle_view_edge_feature_intensities() {       // edge samples
  set_draw_render_types(k_edge_voronoi, k_feature_intensity);
  changed();
}

void Scene::toggle_view_edge_capacities() {
  set_draw_render_types(k_edge_voronoi, k_capacity);
  changed();
}

void Scene::toggle_view_edge_weights() {
  set_draw_render_types(k_edge_voronoi, k_weight);
  changed();
}

void Scene::toggle_view_facet_feature_intensities() {      // facet samples
  set_draw_render_types(k_facet_voronoi, k_feature_intensity);
  changed();
}

void Scene::toggle_view_facet_capacities() {
  set_draw_render_types(k_facet_voronoi, k_capacity);
  changed();
}

void Scene::toggle_view_facet_weights() {
  set_draw_render_types(k_facet_voronoi, k_weight);
  changed();
}

void Scene::eliminate_degenerations() {
  if (m_pInput == NULL) {
    std::cout << "Please open a file first" << std::endl;
    return;
  }
  CGAL::Timer timer;
  timer.start();
  std::cout << std::endl << "Eliminate Input degenerated facets...";
  size_t nb_eliminations = remesh.eliminate_degenerated_facets(m_pInput);
  std::cout << "Done (" << nb_eliminations << " facets eliminated, "
    << timer.time() << " s)" << std::endl;
  if (nb_eliminations > 0) {
    m_pInput->calculate_normals("Input");
    remesh.calculate_feature_intensities("Input", m_pInput);
    m_input_aabb_tree_constructed = false;
    reset();
    reset_draw_render_types();
    changed();
  }
}

void Scene::split_input_long_edges() {
  if (m_pInput == NULL) {
    std::cout << "Please open a file first" << std::endl;
    return;
  }
  CGAL::Timer timer;
  timer.start();
  std::cout << std::endl << "Split long edges for Input...";
  size_t nb_split = remesh.split_long_edges(m_pInput);
  std::cout << "Done (" << nb_split << " edges splited, "
    << timer.time() << " s)" << std::endl;
  if (nb_split > 0) {
    m_pInput->calculate_normals("Input");
    remesh.calculate_feature_intensities("Input", m_pInput);
    m_input_aabb_tree_constructed = false;
    reset();    // do we not need clear links because m_pRemesh is reset
    reset_draw_render_types();
    changed();
  }
}

void Scene::input_properties() {
  if (m_pInput == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_pInput->trace_properties("Input");
  }
}

void Scene::reset_from_input() {
  if (m_pInput == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Resetting from Input...";
    reset();
    std::cout << "Done (" << timer.time() << " s)" << std::endl;
    reset_draw_render_types();
    changed();
  }
}

void Scene::generate_links_and_types() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
    return;
  }
  else {
    generate();
    m_links_initialized = true;
    changed();
  }
}

void Scene::remesh_properties() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      // special usage: we update inside generate_links_and_types()
      generate_links_and_types();
    }
    m_pRemesh->trace_properties("Remesh");  // basic properties
    m_pRemesh->trace_additional_properties(m_pInput->get_diagonal_length());
  }
}

void Scene::isotropic_remeshing() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.isotropic_remeshing(m_input_facet_tree, m_bbox, m_pRemesh);
    changed();
  }
}

void Scene::initial_mesh_simplification() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.initial_mesh_simplification(m_input_facet_tree, m_bbox, m_pRemesh);
    reset_draw_render_types();
    changed();
  }
}

void Scene::split_local_longest_edge() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.split_local_longest_edge(m_input_facet_tree, m_bbox, m_pRemesh);
    reset_draw_render_types();
    changed();
  }
}

void Scene::increase_minimal_angle() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.increase_minimal_angle(m_input_facet_tree, m_bbox, m_pRemesh);
    reset_draw_render_types();
    changed();
  }
}

void Scene::maximize_minimal_angle() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.maximize_minimal_angle(m_input_facet_tree, m_bbox, m_pRemesh);
    reset_draw_render_types();
    changed();
  }
}

void Scene::final_vertex_relocation() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    if (!m_links_initialized) {
      generate();
      m_links_initialized = true;
    }
    std::cout << std::endl;
    remesh.final_vertex_relocation(m_input_facet_tree, m_bbox, m_pRemesh);
    reset_draw_render_types();
    changed();
  }
}

void Scene::update_feature_intensities_and_clear_links() {
  if (m_pInput != NULL && m_pRemesh != NULL) {
    // step 1: update feature intensities
    remesh.calculate_feature_intensities("Input", m_pInput);
    remesh.calculate_feature_intensities("Remesh", m_pRemesh);
    // step 2: clear the links
    remesh.clear_links(m_pInput, m_pRemesh);
    m_links_initialized = false;
    remesh.initialize_private_data();
    changed();
  }
}

void Scene::set_draw_render_types(DrawType draw_type, RenderType render_type) {
  m_draw_type = draw_type;
  m_render_type = render_type;
}

void Scene::reset_draw_render_types() {
  set_draw_render_types(k_polyhedron, k_plain_facets);
}

void Scene::reset() {
  if (m_pInput != NULL) {
    if (m_pRemesh != NULL) {
      delete m_pRemesh;
    }
    m_pRemesh = new Polyhedron(*m_pInput);
    remesh.clear_links(m_pInput, m_pRemesh);
    m_links_initialized = false;
    remesh.initialize_private_data();
  }
}

void Scene::generate() {
  if (m_pInput != NULL && m_pRemesh != NULL) {
    std::cout << std::endl;
    if (!m_input_aabb_tree_constructed) {
      remesh.build_facet_tree(*m_pInput, "Input", &m_input_facet_tree);
      m_input_aabb_tree_constructed = true;
    }
    remesh.build_facet_tree(*m_pRemesh, "Remesh", &m_remesh_facet_tree);
    remesh.generate_links(m_input_facet_tree, m_remesh_facet_tree,
                          m_pInput, m_pRemesh);
  }
}

/*
void Scene::timerEvent(QTimerEvent *) {
  if (manipulatedFrame()->isSpinning())
  set_fast_distance(true);
  ready_to_cut = true;
}

void Scene::draw_vertex_voronoi() {
  bool inherit_element_types = remesh.get_inherit_element_types();
  double sum_theta_value = remesh.get_sum_theta() * CGAL_PI;
  double dihedral_theta_value = remesh.get_dihedral_theta() * CGAL_PI;
  double feature_control_delta = remesh.get_feature_control_delta();
  if (m_view_input && !m_view_remesh) {
    m_pInput->render_vertex_voronoi(m_render_type, inherit_element_types,
      m_view_polyhedron_edges, feature_control_delta, sum_theta_value,
      dihedral_theta_value);
  }
  if (!m_view_input && m_view_remesh) {
    m_pRemesh->render_vertex_voronoi(m_render_type, inherit_element_types,
      m_view_polyhedron_edges, feature_control_delta, sum_theta_value,
      dihedral_theta_value);
  }
}

void Scene::draw_edge_voronoi() {
  double sum_theta_value = remesh.get_sum_theta() * CGAL_PI;
  double dihedral_theta_value = remesh.get_dihedral_theta() * CGAL_PI;
  if (m_view_input && !m_view_remesh) {
    m_pInput->render_edge_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
  if (!m_view_input && m_view_remesh) {
    m_pRemesh->render_edge_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
}

void Scene::draw_facet_voronoi() {
  double sum_theta_value = remesh.get_sum_theta() * CGAL_PI;
  double dihedral_theta_value = remesh.get_dihedral_theta() * CGAL_PI;
  if (m_view_input && !m_view_remesh) {
    m_pInput->render_facet_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
  if (!m_view_input && m_view_remesh) {
    m_pRemesh->render_facet_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
}

void Scene::draw_all_voronoi() {
  double sum_theta_value = remesh.get_sum_theta() * CGAL_PI;
  double dihedral_theta_value = remesh.get_dihedral_theta() * CGAL_PI;
  if (m_view_input && !m_view_remesh) {
    m_pInput->render_all_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
  if (!m_view_input && m_view_remesh) {
    m_pRemesh->render_all_voronoi(m_render_type, m_view_polyhedron_edges,
      sum_theta_value, dihedral_theta_value);
  }
}*/