#include "Scene.h"

//#include <fstream>
//#include <QFileInfo>
//#include <QOpenGLShader>
//#include <QDebug>
//#include <CGAL/Polygon_mesh_processing/border.h>
//#include <CGAL/IO/PLY_reader.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <boost/algorithm/string.hpp>
// TODO(kaimo hu): replace line 13 with line 12 when published
// #include <CGAL/Polygon_mesh_processing/minimal_angle_remeshing.h>
#include "../../include/CGAL/Polygon_mesh_processing/minimal_angle_remeshing.h"

Scene::Scene() 
  : m_frame(new ManipulatedFrame()), 
    gl_init(false),
    are_buffers_initialized(false),
    gl(NULL),
    m_pInput(NULL),
    m_pRemesh(NULL) {
  // 1) member data initialization
  startTimer(0);
  m_bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX, 
                DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
  // 2) view option initialization
  m_view_input = false;
  m_view_remesh = true;
  m_view_face_in_start_points = false;
  m_view_face_in_end_points = false;
  m_view_face_in_links = false;
  m_view_face_out_start_points= false;
  m_view_face_out_end_points = false;
  m_view_face_out_links = false;
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
  m_view_mesh_edges = true;
  reset_draw_render_types();
  // 3) parameters initialization
  m_target_edge_length = TARGET_EDGE_LENGTH;
  m_smooth_iteration_count = SMOOTH_ITERATION_COUNT;
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
  std::cout << "Computing bbox...";
  m_bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
                DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
  if (m_pInput != NULL && !m_pInput->is_empty()) {
    const Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
    Bbox bbox_input = remesher->get_input_bbox();
    m_bbox = m_bbox + bbox_input;
    std::cout << "Done" << std::endl;
  }
  else {
    std::cout << "failed (no surface mesh or empty surface mesh)." << std::endl;
  }
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
  // input faces
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

  // remesh faces
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
  vbo[VBOs::kb_Remesh_face_colors].release();
  vao[VAOs::ka_Remesh_faces].release();

  // Input cell boundaries
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

  // Remesh cell boundaries
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

  // Face in start points
  vao[VAOs::ka_Face_in_start].bind();
  vbo[VBOs::kb_Face_in_start_points].bind();
  vbo[VBOs::kb_Face_in_start_points].allocate(
      pos_face_in_start_points.data(),
      static_cast<int>(pos_face_in_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Face_in_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Face_in_start].release();

  // Face in end points
  vao[VAOs::ka_Face_in_end].bind();
  vbo[VBOs::kb_Face_in_end_points].bind();
  vbo[VBOs::kb_Face_in_end_points].allocate(pos_face_in_end_points.data(),
      static_cast<int>(pos_face_in_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Face_in_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Face_in_end].release();

  // Face in links
  vao[VAOs::ka_Face_in_links].bind();
  vbo[VBOs::kb_Face_in_link_lines].bind();
  vbo[VBOs::kb_Face_in_link_lines].allocate(pos_face_in_links.data(),
      static_cast<int>(pos_face_in_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Face_in_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Face_in_links].release();

  // Face out start points
  vao[VAOs::ka_Face_out_start].bind();
  vbo[VBOs::kb_Face_out_start_points].bind();
  vbo[VBOs::kb_Face_out_start_points].allocate(
      pos_face_out_start_points.data(),
      static_cast<int>(pos_face_out_start_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Face_out_start_points].release();
  rendering_program.release();
  vao[VAOs::ka_Face_out_start].release();

  // Face out end points
  vao[VAOs::ka_Face_out_end].bind();
  vbo[VBOs::kb_Face_out_end_points].bind();
  vbo[VBOs::kb_Face_out_end_points].allocate(pos_face_out_end_points.data(),
      static_cast<int>(pos_face_out_end_points.size() * sizeof(float)));
  points_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(points_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(points_vertexLocation);
  vbo[VBOs::kb_Face_out_end_points].release();
  rendering_program.release();
  vao[VAOs::ka_Face_out_end].release();

  // Face out links
  vao[VAOs::ka_Face_out_links].bind();
  vbo[VBOs::kb_Face_out_link_lines].bind();
  vbo[VBOs::kb_Face_out_link_lines].allocate(pos_face_out_links.data(),
      static_cast<int>(pos_face_out_links.size() * sizeof(float)));
  lines_vertexLocation = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.setAttributeBuffer(lines_vertexLocation, GL_FLOAT, 0, 3);
  rendering_program.enableAttributeArray(lines_vertexLocation);
  vbo[VBOs::kb_Face_out_link_lines].release();
  rendering_program.release();
  vao[VAOs::ka_Face_out_links].release();

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
      (m_render_type == k_classifications || m_view_mesh_edges)) {
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
      (m_view_minimal_angle && m_render_type == k_plain_faces))) {
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
      (m_render_type == k_classifications || m_view_mesh_edges)) {
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
      (m_view_minimal_angle && m_render_type == k_plain_faces))) {
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
  if (m_draw_type == DrawType::k_mesh) {
    gl->glDisable(GL_LIGHTING);
    gl->glPointSize(3.0f);
    gl->glLineWidth(1.0f);
    // face in start points
    if (m_view_face_in_start_points && pos_face_in_start_points.size() > 0) {
      vao[VAOs::ka_Face_in_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 250, 250);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_face_in_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_in_start].release();
    }
    // face in links
    if (m_view_face_in_links && pos_face_in_links.size() > 0) {
      vao[VAOs::ka_Face_in_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 200, 200);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_face_in_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_in_links].release();
    }
    // face in end points
    if (m_view_face_in_end_points && pos_face_in_end_points.size() > 0) {
      vao[VAOs::ka_Face_in_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(0, 150, 150);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_face_in_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_in_end].release();
    }
    // face out start points
    if (m_view_face_out_start_points && 
        pos_face_out_start_points.size() > 0) {
      vao[VAOs::ka_Face_out_start].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(250, 0, 250);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_face_out_start_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_out_start].release();
    }
    // face out links
    if (m_view_face_out_links && pos_face_out_links.size() > 0) {
      vao[VAOs::ka_Face_out_links].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(200, 0, 200);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_LINES, 0,
        static_cast<GLsizei>(pos_face_out_links.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_out_links].release();
    }
    // face out end points
    if (m_view_face_out_end_points && pos_face_out_end_points.size() > 0) {
      vao[VAOs::ka_Face_out_end].bind();
      attrib_buffers(viewer);
      rendering_program.bind();
      color.setRgb(150, 0, 150);
      rendering_program.setUniformValue(colorLocation, color);
      rendering_program.setUniformValue(fLocation, fMatrix);
      gl->glDrawArrays(GL_POINTS, 0,
        static_cast<GLsizei>(pos_face_out_end_points.size() / 3));
      rendering_program.release();
      vao[VAOs::ka_Face_out_end].release();
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
    // step 1: compute the faces
    Color face_color(150, 150, 200);
    compute_faces(true, face_color, &pos_input_faces, &pos_input_face_normals,
        &pos_input_face_colors, &pos_input_boundaries, &pos_input_samples);
    // step 2: compute the edges
    if (m_render_type == k_classifications) {
      compute_classified_edges(true, &pos_input_normal_edges,
                               &pos_input_special_edges);
    }
    else {
      if (m_view_mesh_edges) {
        compute_edges(true, &pos_input_normal_edges);
      }
      if (m_view_minimal_angle && m_render_type == k_plain_faces) {
        compute_min_radian_edges(true, &pos_input_special_edges);
      }
    }
  }
  if (m_view_remesh && m_pRemesh != NULL) {
    // step 1: compute the faces
    Color face_color(200, 150, 150);
    compute_faces(false, face_color, &pos_remesh_faces, 
        &pos_remesh_face_normals, &pos_remesh_face_colors,
        &pos_remesh_boundaries, &pos_remesh_samples);
    // step 2: compute the edges
    if (m_render_type == k_classifications) {
      compute_classified_edges(false, &pos_remesh_normal_edges,
        &pos_remesh_special_edges);
    }
    else {
      if (m_view_mesh_edges) {
        compute_edges(false, &pos_remesh_normal_edges);
      }
      if (m_view_minimal_angle && m_render_type == k_plain_faces) {
        compute_min_radian_edges(false, &pos_remesh_special_edges);
      }
    }
  }
  // step 3: compute samples and links
  if (m_pInput != NULL && m_pRemesh != NULL
    && m_draw_type == DrawType::k_mesh) {
    // input
    if (m_view_face_in_start_points) {    // face in links
      compute_face_start_points(true, &pos_face_in_start_points);
    }
    if (m_view_face_in_end_points) {
      compute_face_end_points(true, &pos_face_in_end_points);
    }
    if (m_view_face_in_links) {
      compute_face_links(true, &pos_face_in_links);
    }
    if (m_view_edge_in_start_points) {    // edge in links
      compute_edge_start_points(true, &pos_edge_in_start_points);
    }
    if (m_view_edge_in_end_points) {
      compute_edge_end_points(true, &pos_edge_in_end_points);
    }
    if (m_view_edge_in_links) {
      compute_edge_links(true, &pos_edge_in_links);
    }
    if (m_view_vertex_in_start_points) {  // vertex in links
      compute_vertex_start_points(true, &pos_vertex_in_start_points);
    }
    if (m_view_vertex_in_end_points) {
      compute_vertex_end_points(true, &pos_vertex_in_end_points);
    }
    if (m_view_vertex_in_links) {
      compute_vertex_links(true, &pos_vertex_in_links);
    }
    // remesh
    if (m_view_face_out_start_points) {   // face out links
      compute_face_start_points(false, &pos_face_out_start_points);
    }
    if (m_view_face_out_end_points) {
      compute_face_end_points(false, &pos_face_out_end_points);
    }
    if (m_view_face_out_links) {
      compute_face_links(false, &pos_face_out_links);
    }
    if (m_view_edge_out_start_points) {   // edge out links
      compute_edge_start_points(false, &pos_edge_out_start_points);
    }
    if (m_view_edge_out_end_points) {
      compute_edge_end_points(false, &pos_edge_out_end_points);
    }
    if (m_view_edge_out_links) {
      compute_edge_links(false, &pos_edge_out_links);
    }
    if (m_view_vertex_out_start_points) { // vertex out links
      compute_vertex_start_points(false, &pos_vertex_out_start_points);
    }
    if (m_view_vertex_out_end_points) {
      compute_vertex_end_points(false, &pos_vertex_out_end_points);
    }
    if (m_view_vertex_out_links) {
      compute_vertex_links(false, &pos_vertex_out_links);
    }
  }
}

void Scene::compute_faces(bool is_input, Color face_color, 
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
    std::vector<float> *pos_samples) const {
  FT sum_theta_value = get_sum_theta() * CGAL_PI;
  FT dihedral_theta_value = get_dihedral_theta() * CGAL_PI;
  const Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
  FT max_error_threshold_value = remesher->get_max_error_threshold_value();
  pos_faces->resize(0);
  pos_face_normals->resize(0);
  pos_face_colors->resize(0);
  pos_boundaries->resize(0);
  pos_samples->resize(0);
  switch (m_draw_type) {
  case DrawType::k_mesh:
    compute_mesh_faces(is_input, face_color, sum_theta_value,
        dihedral_theta_value, max_error_threshold_value, pos_faces,
        pos_face_normals, pos_face_colors);
    break;
  case DrawType::k_all_voronoi:
    compute_all_voronois(is_input, sum_theta_value, dihedral_theta_value,
        pos_faces, pos_face_normals, pos_face_colors, pos_boundaries,
        pos_samples);
    break;
  case DrawType::k_vertex_voronoi:
    compute_vertex_voronois(is_input, sum_theta_value, dihedral_theta_value, 
        pos_faces, pos_face_normals, pos_face_colors, pos_boundaries,
        pos_samples);
    break;
  case DrawType::k_edge_voronoi:
    compute_edge_voronois(is_input, sum_theta_value, dihedral_theta_value,
        pos_faces, pos_face_normals, pos_face_colors, pos_boundaries,
        pos_samples);
    break;
  case DrawType::k_face_voronoi:
    compute_face_voronois(is_input, sum_theta_value, dihedral_theta_value,
        pos_faces, pos_face_normals, pos_face_colors, pos_boundaries,
        pos_samples);
    break;
  default:
    break;
  }
}

void Scene::compute_mesh_faces(bool is_input, Color face_color,
    FT sum_theta_value, FT dihedral_theta_value, FT max_error_threshold_value,
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors) const {
  FT min_value = 0;
  FT max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    face_descriptor fd = *fi;
    const Normal &normal = mesh_properties->get_face_normal(fd);
    halfedge_descriptor hd = mesh.halfedge(fd);
    vertex_descriptor vd = mesh.target(hd);
    FT error = CGAL::sqrt(mesh_properties->get_face_max_squared_error(fd));
    if (error > max_error_threshold_value) {
      error = max_error_threshold_value;
    }
    else {
      // if error == max_error_threshold_value, normalize it to 0.5
      error = error / (2 * max_error_threshold_value);  // normalization
    }
    Color error_color = get_rgb_color(0, error, 1.0);
    Color ifi_color;
    for (int i = 0; i <= 2; ++i) {
      compute_vertex(mesh, vd, pos_faces);
      compute_normal(normal, pos_face_normals);
      switch (m_render_type) {
      case RenderType::k_plain_faces:
        compute_color(face_color, pos_face_colors);
        break;
      case RenderType::k_ifi_faces:
        ifi_color = get_vertex_sample_normalized_color(is_input, vd, min_value,
                                                       max_value, 240);
        compute_color(ifi_color, pos_face_colors);
        break;
      case RenderType::k_mr_faces:
        compute_color(error_color, pos_face_colors);
        break;
      default:
        break;
      }
      hd = mesh.next(hd);
      vd = mesh.target(hd);
    }
  }
}

void Scene::compute_all_voronois(bool is_input, FT sum_theta_value,
    FT dihedral_theta_value, std::vector<float> *pos_faces,
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries, 
    std::vector<float> *pos_samples) const {
  // step 1: compute all the minimal value and maximal value
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  FT min_value = 0.0, max_value = 0.0;
  switch (m_render_type) {
  case RenderType::k_feature_intensity:
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
    break;
  case RenderType::k_capacity:
    max_value = 2.0;
    break;
  case RenderType::k_weight:
    max_value = mesh_properties->get_max_sample_weight();
    break;
  default:
    break;
  }
  // step 2: compute all sample cells and cell boundaries
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    const Normal &normal = mesh_properties->get_face_normal(*fi);
    Point_list samples;
    Color_list colors;
    get_all_sample_normalized_colors(is_input, *fi, min_value, max_value,
                                     240, &samples, &colors);
    Bvd bvd(mesh_properties->triangle(*fi));
    bvd.compute_voronoi_cells_and_boundaries(samples, normal, colors,
      pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
  }
  // step 3: compute all samples
  compute_all_samples(is_input, pos_samples);
}

void Scene::compute_vertex_voronois(bool is_input, FT sum_theta_value, 
    FT dihedral_theta_value, std::vector<float> *pos_faces,
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries,
    std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  // step 1: comput the minimal and maximal value
  FT min_value = 0.0, max_value = 0.0;
  switch (m_render_type) {
  case RenderType::k_classifications:
    break;
  case RenderType::k_gaussian_curvature:
    max_value = sum_theta_value;
    break;
  case RenderType::k_maximal_halfedge_dihedral:
    max_value = dihedral_theta_value;
    break;
  case RenderType::k_feature_intensity:
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
    break;
  case RenderType::k_capacity:
    //max_value = 2 * get_max_sample_capacity();
    max_value = 2 * mesh_properties->get_max_vertex_sample_capacity();
    break;
  case RenderType::k_weight:
    //max_value = get_max_sample_weight();
    max_value = mesh_properties->get_max_vertex_sample_weight();
    break;
  default:
    break;
  }
  // step 2: compute the cell and cell boundaries
  Color color;
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
    vi != mesh.vertices().end(); ++vi) {
    if (m_render_type == RenderType::k_classifications) {
      color = get_vertex_classification_color(is_input, *vi);
    }
    else {
      color = get_vertex_sample_normalized_color(is_input, *vi, min_value,
                                                 max_value, 240);
    }
    const Point &p = mesh.point(*vi);
    Halfedge_around_target_circulator hb(mesh.halfedge(*vi), mesh), he(hb);
    do {
      if (!mesh.is_border(*hb)) {
        // step 1: compute the face
        face_descriptor fd = mesh.face(*hb);
        const Normal &normal = mesh_properties->get_face_normal(fd);
        const Point &p1 = mesh_properties->midpoint(*hb);
        const Point c = mesh_properties->centroid(fd);
        const Point &p2 = mesh_properties->midpoint(mesh.next(*hb));
        compute_triangle(p1, p, c, normal, color,
          pos_faces, pos_face_normals, pos_face_colors);
        compute_triangle(c, p, p2, normal, color,
          pos_faces, pos_face_normals, pos_face_colors);
        // step 2: compute the edges
        compute_segment(p1, c, pos_boundaries);
        compute_segment(c, p2, pos_boundaries);
      }
      ++hb;
    } while (hb != he);
  }
  // step 3: compute vertices
  compute_vertices(is_input, pos_samples);
}

void Scene::compute_edge_voronois(bool is_input, FT sum_theta_value,
    FT dihedral_theta_value, std::vector<float> *pos_faces,
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries,
    std::vector<float> *pos_samples) const {
  if (m_render_type == RenderType::k_normal_dihedral) {
    // only render two triangles for each edge, Do not render samples
    compute_edge_normal_dihedrals(is_input, dihedral_theta_value, pos_faces,
        pos_face_normals, pos_face_colors, pos_boundaries);
  }
  else {
    // we need to render triangles for each sample, render the samples
    compute_edge_sample_properties(is_input, sum_theta_value,
        dihedral_theta_value, pos_faces, pos_face_normals, pos_face_colors,
        pos_boundaries, pos_samples);
  }
}

void Scene::compute_face_voronois(bool is_input, FT sum_theta_value,
    FT dihedral_theta_value, std::vector<float> *pos_faces,
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries, 
    std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  // step 1: compute the minimal and maximal values
  FT min_value = 0.0, max_value = 0.0;
  switch (m_render_type) {
  case RenderType::k_feature_intensity:
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
    break;
  case RenderType::k_capacity:
    //max_value = 2 * get_max_sample_capacity();
    max_value = 2 * mesh_properties->get_max_face_sample_capacity();
    break;
  case RenderType::k_weight:
    //max_value = get_max_sample_weight();
    max_value = mesh_properties->get_max_face_sample_weight();
    break;
  default:
    break;
  }
  // step 2: compute face_voronoi cell and boundaries
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    const Normal &normal = mesh_properties->get_face_normal(*fi);
    Point_list samples;
    Color_list colors;
    get_face_sample_normalized_colors(is_input, *fi, min_value, 
                                      max_value, 240, &samples, &colors);
    Bvd bvd(mesh_properties->triangle(*fi));
    bvd.compute_voronoi_cells_and_boundaries(samples, normal, colors,
      pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
  }
  // step 3: compute samples
  compute_face_samples(is_input, pos_samples);
}

void Scene::compute_edge_normal_dihedrals(bool is_input,
    FT dihedral_theta_value, std::vector<float> *pos_faces, 
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  FT min_value = 0.0, max_value = dihedral_theta_value;
  Color color;
  // step 1: compute the faces
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
      ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    FT value = mesh_properties->get_halfedge_normal_dihedral(hd) / max_value;
    color = get_rgb_color(240, value, 1.0);
    face_descriptor fd = mesh.face(hd);            // the first triangle
    const Normal &normal = mesh_properties->get_face_normal(fd);
    const Point &a = mesh.point(mesh.target(hd));
    const Point b = mesh_properties->centroid(fd);
    const Point &c = mesh.point(mesh.source(hd));
    compute_triangle(a, b, c, normal, color,
      pos_faces, pos_face_normals, pos_face_colors);
    if (!mesh.is_border(mesh.opposite(hd))) {     // the second triangle
      hd = mesh.opposite(hd);
      fd = mesh.face(hd);
      const Normal &normal = mesh_properties->get_face_normal(fd);
      const Point &a = mesh.point(mesh.target(hd));
      const Point b = mesh_properties->centroid(fd);
      const Point &c = mesh.point(mesh.source(hd));
      compute_triangle(a, b, c, normal, color, pos_faces, pos_face_normals, 
                       pos_face_colors);
    }
  }
  // step 2: compute the edges
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
    fi != mesh.faces().end(); ++fi) {
    halfedge_descriptor hd = mesh.halfedge(*fi);
    const Point &a = mesh.point(mesh.source(hd));
    const Point &b = mesh.point(mesh.target(hd));
    const Point &c = mesh.point(mesh_properties->get_opposite_vertex(hd));
    Point d = mesh_properties->centroid(*fi);
    compute_segment(d, a, pos_boundaries);
    compute_segment(d, b, pos_boundaries);
    compute_segment(d, c, pos_boundaries);
  }
}

void Scene::compute_edge_sample_properties(bool is_input, FT sum_theta_value,
    FT dihedral_theta_value, std::vector<float> *pos_faces, 
    std::vector<float> *pos_face_normals, std::vector<float> *pos_face_colors,
    std::vector<float> *pos_boundaries,
    std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  // step 1: get the minimal and maximal values
  FT min_value = 0.0, max_value = 0.0;
  switch (m_render_type) {
  case RenderType::k_feature_intensity:
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
    break;
  case RenderType::k_capacity:
    //max_value = 2 * get_max_sample_capacity();
    max_value = 2 * mesh_properties->get_max_edge_sample_capacity();
    break;
  case RenderType::k_weight:
    //max_value = get_max_sample_weight();
    max_value = mesh_properties->get_max_edge_sample_weight();
    break;
  default:
    break;
  }
  // step 2: compute the faces and faces
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    // step 2.1: compute the samples
    const Link_list &edge_out_links = 
        mesh_properties->get_halfedge_out_links(hd);
    for (auto cit = edge_out_links.begin();
      cit != edge_out_links.end(); ++cit) {
      const Link &link = *cit;
      compute_point(link.second.first, pos_samples);
    }
    // step 2.2: compute the faces and edges
    Link_list_const_iter first = edge_out_links.cbegin(), second = first;
    ++second;
    Point_list points;
    points.push_back(mesh.point(mesh.source(hd)));
    for (; second != edge_out_links.cend(); ++first, ++second) {
      points.push_back(CGAL::midpoint(first->second.first,
        second->second.first));
    }
    points.push_back(mesh.point(mesh.target(hd)));
    Color_list colors;
    get_edge_sample_normalized_colors(is_input, hd, min_value, 
                                      max_value, 240, &colors);
    face_descriptor fd = mesh.face(hd);        // first triangle
    Point c = mesh_properties->centroid(fd);
    const Normal &normal = mesh_properties->get_face_normal(fd);
    Point_iter it1 = points.begin(), it2 = it1;
    ++it2;
    Color_const_iter cit = colors.begin();
    for (; it2 != points.end(); ++it1, ++it2, ++cit) {
      if (it1 == points.begin()) {
        compute_segment(c, *it1, pos_boundaries);
      }
      compute_triangle(*it1, *it2, c, normal, *cit, pos_faces,
                       pos_face_normals, pos_face_colors);
    }
    if (!mesh.is_border(mesh.opposite(hd))) {  // second triangle
      fd = mesh.face(mesh.opposite(hd));
      c = mesh_properties->centroid(fd);
      const Normal &normal = mesh_properties->get_face_normal(fd);
      it1 = points.begin(), it2 = it1;
      ++it2;
      cit = colors.begin();
      for (; it2 != points.end(); ++it1, ++it2, ++cit) {
        if (it1 != points.begin()) {
          compute_segment(c, *it1, pos_boundaries);
        }
        compute_triangle(*it1, *it2, c, normal, *cit,
          pos_faces, pos_face_normals, pos_face_colors);
      }
    }
  }
  // step 3: compute the rest of the edges
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    Point c = mesh_properties->centroid(*fi);
    halfedge_descriptor hd = mesh.halfedge(*fi);
    for (int i = 0; i <= 2; ++i) {
      compute_segment(c, mesh.point(mesh.target(hd)), pos_boundaries);
      hd = mesh.next(hd);
    }
  }
}

Color Scene::get_vertex_classification_color(bool is_input,
                                             vertex_descriptor vd) const {
  Halfedge_list effective_edges;
  const Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  VertexType vertex_type = mesh_properties->get_vertex_type(
      vd, &effective_edges, remesher->get_named_parameters());
  if (vertex_type == VertexType::k_feature_vertex) {
    return get_rgb_color(300, 0.8, 1.0);      // feature vertex: purple
  }
  else if (vertex_type == VertexType::k_crease_vertex) {
    return get_rgb_color(180, 0.8, 1.0);      // crease vertex: cyan
  }
  else {
    return get_rgb_color(60, 0.8, 1.0);       // smooth vertex: yellow
  }
}

void Scene::compute_all_samples(bool is_input,
                                std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  // step 1: face out samples
  compute_face_samples(is_input, pos_samples);
  // step 2: edge out samples
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    const Link_list &edge_out_links =
        mesh_properties->get_halfedge_out_links(hd);
    for (auto it = edge_out_links.begin();
      it != edge_out_links.end(); ++it) {
      const Point &p = it->second.first;
      compute_point(p, pos_samples);
    }
  }
  // step 3: vertex out samples
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
    vi != mesh.vertices().end(); ++vi) {
    const Link &link = mesh_properties->get_vertex_out_link(*vi);
    const Point &p = link.second.first;
    compute_point(p, pos_samples);
  }
}

void Scene::compute_face_samples(bool is_input,
                                 std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
    fi != mesh.faces().end(); ++fi) {
    const Link_list &face_out_links = mesh_properties->get_face_out_links(*fi);
    for (auto it = face_out_links.begin(); it != face_out_links.end(); ++it) {
      const Link &link = *it;
      compute_point(link.second.first, pos_samples);
    }
  }
}

Color Scene::get_vertex_sample_normalized_color(bool is_input,
    vertex_descriptor vd, FT min_value, FT max_value, FT h) const {
  FT value = 0.0;
  FT weight = 0.0, feature_intensity = 0.0;
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  switch (m_render_type) {
  case RenderType::k_gaussian_curvature:
    value = mesh_properties->get_vertex_gaussian_curvature(vd);
    break;
  case RenderType::k_maximal_halfedge_dihedral:
    value = mesh_properties->get_vertex_max_dihedral(vd);
    break;
  case RenderType::k_ifi_faces:
  case RenderType::k_feature_intensity:
    value = mesh_properties->calculate_feature_intensity(vd) - 1.0;
    break;
  case RenderType::k_capacity:
    weight = mesh_properties->get_vertex_out_link(vd).first;
    feature_intensity = mesh_properties->calculate_feature_intensity(vd);
    value = weight / feature_intensity;
    break;
  case RenderType::k_weight:
    value = mesh_properties->get_vertex_out_link(vd).first;
    break;
  default:
    break;
  }
  if (value > max_value) {
    value = max_value;
  }
  if (value < min_value) {
    value = min_value;
  }
  value = (value - min_value) / (max_value - min_value);  // normalize
  return get_rgb_color(h, value, 1.0);
}

void Scene::get_all_sample_normalized_colors(bool is_input, 
    face_descriptor fd, FT min_value, FT max_value, FT h, Point_list *samples,
    Color_list *colors) const {
  // calculate the samples and their weights
  std::map<Point, double, Point_Comp> disturbed_border_samples;
  get_disturbed_border_samples_with_weights(is_input, fd, 0.01,
                                            &disturbed_border_samples);
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Link_list &face_out_links = mesh_properties->get_face_out_links(fd);
  for (auto it = face_out_links.begin(); it != face_out_links.end(); ++it) {
    const Link &link = *it;
    if (disturbed_border_samples.find(link.second.first) ==
      disturbed_border_samples.end()) {
      switch (m_render_type) {
      case RenderType::k_feature_intensity:
        disturbed_border_samples[link.second.first] = link.first - 1.0;
        break;
      case RenderType::k_weight:
        disturbed_border_samples[link.second.first] = link.first;
        break;
      case RenderType::k_capacity:
        disturbed_border_samples[link.second.first] = 1.0;
        break;
      default:
        break;
      }
    }
  }
  // normalize the weights, and convert to colors
  for (auto it = disturbed_border_samples.begin();
    it != disturbed_border_samples.end(); ++it) {
    if (it->second > max_value) {
      it->second = max_value;
    }
    if (it->second < min_value) {
      it->second = min_value;
    }
    it->second = (it->second - min_value) / (max_value - min_value);
    samples->push_back(it->first);
    colors->push_back(get_rgb_color(h, it->second, 1.0));
  }
}

void Scene::get_edge_sample_normalized_colors(bool is_input,
    halfedge_descriptor hd, FT min_value, FT max_value, FT h,
    Color_list *colors) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  FT face_area = mesh_properties->area(mesh.face(hd));   // the face area
  if (!mesh.is_border(mesh.opposite(hd))) {
    face_area += mesh_properties->area(mesh.face(mesh.opposite(hd)));
  }
  const Link_list &edge_out_links =
      mesh_properties->get_halfedge_out_links(hd);
  FT capacity = face_area / (3 * edge_out_links.size());
  std::list<double> values;
  // step 1: calculate the original values
  for (auto it = edge_out_links.begin(); it != edge_out_links.end(); ++it) {
    switch (m_render_type) {
    case RenderType::k_feature_intensity:
      values.push_back(it->first / capacity - 1.0);
      break;
    case RenderType::k_capacity:
      values.push_back(capacity);
      break;
    case RenderType::k_weight:
      values.push_back(it->first);
      break;
    default:
      break;
    }
  }
  // step 2: normalize the values, and convert to colors
  for (auto it = values.begin(); it != values.end(); ++it) {
    if (*it > max_value) {
      *it = max_value;
    }
    if (*it < min_value) {
      *it = min_value;
    }
    *it = (*it - min_value) / (max_value - min_value);
    colors->push_back(get_rgb_color(h, *it, 1.0));
  }
}

void Scene::get_face_sample_normalized_colors(bool is_input,
    face_descriptor fd, FT min_value, FT max_value, FT h, Point_list *samples,
    Color_list *colors) const {
  // step 1: calculate the original values
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  std::list<double> values;
  const Link_list &face_out_links = mesh_properties->get_face_out_links(fd);
  FT capacity = mesh_properties->area(fd) / face_out_links.size();
  for (auto it = face_out_links.begin(); it != face_out_links.end(); ++it) {
    samples->push_back(it->second.first);
    switch (m_render_type) {
    case RenderType::k_feature_intensity:
      values.push_back(it->first / capacity - 1.0);
      break;
    case RenderType::k_capacity:
      values.push_back(capacity);
      break;
    case RenderType::k_weight:
      values.push_back(it->first);
      break;
    default:
      break;
    }
  }
  // step 2: normalize the values, and convert to colors
  for (auto it = values.begin(); it != values.end(); ++it) {
    if (*it > max_value) {
      *it = max_value;
    }
    if (*it < min_value) {
      *it = min_value;
    }
    *it = (*it - min_value) / (max_value - min_value);
    colors->push_back(get_rgb_color(h, *it, 1.0));
  }
}

void Scene::get_disturbed_border_samples_with_weights(bool is_input,
    face_descriptor fd, FT disturb_ratio,
    std::map<Point, double, Point_Comp> *disturbed_border_samples) const {
  // this version is used for rendering
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  halfedge_descriptor hd = mesh.halfedge(fd);
  FT value = 0.0;
  for (int i = 0; i <= 2; ++i) {
    // step 1: add the vertex sample
    vertex_descriptor vd = mesh_properties->get_opposite_vertex(hd);
    Vector vec = mesh_properties->midpoint(hd) - mesh.point(vd);
    vec = vec * disturb_ratio;
    const Link &link = mesh_properties->get_vertex_out_link(vd);
    if (m_render_type == RenderType::k_capacity) {
      value = 1.0;
    }
    else if (m_render_type == RenderType::k_feature_intensity) {
      value = link.first - 1.0;
    }
    else {
      value = link.first;
    }
    disturbed_border_samples->insert(std::make_pair(link.second.first + vec,
                                     value));
    // step 2: add the edge samples
    halfedge_descriptor h = hd;
    if (mesh_properties->get_halfedge_normal_dihedral(h) == -1.0) {
      h = mesh.opposite(h);
    }
    const Link_list &edge_out_links = 
        mesh_properties->get_halfedge_out_links(h);
    for (Link_list_const_iter cit = edge_out_links.begin();
      cit != edge_out_links.end(); ++cit) {
      const Point &p = cit->second.first;
      vec = mesh.point(vd) - p;
      vec = vec * disturb_ratio;
      if (m_render_type == RenderType::k_capacity) {
        value = 1.0;
      }
      else if (m_render_type == RenderType::k_feature_intensity) {
        value = cit->first - 1.0;
      }
      else {
        value = cit->first;
      }
      disturbed_border_samples->insert(std::make_pair(p + vec, value));
    }
    hd = mesh.next(hd);
  }
}

void Scene::compute_classified_edges(bool is_input,
    std::vector<float> *pos_normal_edges,
    std::vector<float> *pos_special_edges) const {
  pos_normal_edges->resize(0);
  pos_special_edges->resize(0);
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    if (mesh_properties->get_halfedge_is_crease(hd)) {
      compute_halfedge(mesh, hd, pos_special_edges);
    }
    else {
      compute_halfedge(mesh, hd, pos_normal_edges);
    }
  }
}

void Scene::compute_vertices(bool is_input,
                             std::vector<float> *pos_samples) const {
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
    vi != mesh.vertices().end(); ++vi) {
    const Point &p = mesh.point(*vi);
    compute_point(p, pos_samples);
  }
}

void Scene::compute_edges(bool is_input, std::vector<float> *pos_edges) const {
  pos_edges->resize(0);
  const Mesh_properties* mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
      ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    compute_halfedge(mesh, hd, pos_edges);
  }
}

void Scene::compute_min_radian_edges(bool is_input,
    std::vector<float> *pos_min_radian_edges) const {
  pos_min_radian_edges->resize(0);
  FT min_radian = CGAL_PI;
  const Mesh_properties* mesh_properties = get_mesh_properties(is_input);
  halfedge_descriptor min_radian_hd = 
      mesh_properties->calculate_minimal_radian(&min_radian);
  if (min_radian != CGAL_PI) {
    const Mesh &mesh = mesh_properties->get_mesh();
    compute_halfedge(mesh, mesh.next(min_radian_hd), pos_min_radian_edges);
    compute_halfedge(mesh, mesh.prev(min_radian_hd), pos_min_radian_edges);
  }
}

void Scene::compute_halfedge(const Mesh &mesh, halfedge_descriptor hd,
                                    std::vector<float> *pos) const {
  const Point &p = mesh.point(mesh.source(hd));
  const Point &q = mesh.point(mesh.target(hd));
  compute_segment(p, q, pos);
}

void Scene::compute_vertex(const Mesh &mesh, vertex_descriptor vd,
                                  std::vector<float> *pos) const {
  compute_point(mesh.point(vd), pos);
}

void Scene::compute_triangle(const Point &p1, const Point &p2, const Point &p3,
    const Normal &normal, const Color &color, std::vector<float> *pos_faces,
    std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors) const {
  compute_triangle_point(p1, normal, color,
    pos_faces, pos_face_normals, pos_face_colors);
  compute_triangle_point(p2, normal, color,
    pos_faces, pos_face_normals, pos_face_colors);
  compute_triangle_point(p3, normal, color,
    pos_faces, pos_face_normals, pos_face_colors);
}

void Scene::compute_triangle_point(const Point &p, const Normal &normal,
  const Color &color, std::vector<float> *pos_faces,
  std::vector<float> *pos_face_normals,
  std::vector<float> *pos_face_colors) const {
  compute_point(p, pos_faces);
  compute_normal(normal, pos_face_normals);
  compute_color(color, pos_face_colors);
}

void Scene::compute_face_start_points(bool is_input,
    std::vector<float> *pos_face_start_point) const {
  pos_face_start_point->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    const Link_list &face_out_links = mesh_properties->get_face_out_links(*fi);
    for (auto it = face_out_links.begin();
      it != face_out_links.end(); ++it) {
      const Link &link = *it;
      compute_point(link.second.first, pos_face_start_point);
    }
  }
}

void Scene::compute_face_end_points(bool is_input,
    std::vector<float> *pos_face_end_point) const {
  pos_face_end_point->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
      fi != mesh.faces().end(); ++fi) {
    const Link_list &face_out_links = mesh_properties->get_face_out_links(*fi);
    for (auto it = face_out_links.begin();
      it != face_out_links.end(); ++it) {
      const Link &link = *it;
      compute_point(link.second.second, pos_face_end_point);
    }
  }
}

void Scene::compute_face_links(bool is_input,
    std::vector<float> *pos_face_links) const {
  pos_face_links->clear();
  const Mesh_properties* mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Face_range::const_iterator fi = mesh.faces().begin();
    fi != mesh.faces().end(); ++fi) {
    const Link_list &face_out_links = mesh_properties->get_face_out_links(*fi);
    for (auto it = face_out_links.begin();
      it != face_out_links.end(); ++it) {
      const Link &link = *it;
      const Point &p = link.second.first;
      const Point &q = link.second.second;
      compute_segment(p, q, pos_face_links);
    }
  }
}

void Scene::compute_edge_start_points(bool is_input,
    std::vector<float> *pos_edge_start_points) const {
  pos_edge_start_points->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    const Link_list &edge_out_links = 
        mesh_properties->get_halfedge_out_links(hd);
    for (auto it = edge_out_links.begin();
      it != edge_out_links.end(); ++it) {
      const Link &link = *it;
      compute_point(link.second.first, pos_edge_start_points);
    }
  }
}

void Scene::compute_edge_end_points(bool is_input,
    std::vector<float> *pos_edge_end_points) const {
  pos_edge_end_points->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    const Link_list &edge_out_links = 
        mesh_properties->get_halfedge_out_links(hd);
    for (auto it = edge_out_links.begin();
      it != edge_out_links.end(); ++it) {
      const Link &link = *it;
      compute_point(link.second.second, pos_edge_end_points);
    }
  }
}

void Scene::compute_edge_links(bool is_input,
    std::vector<float> *pos_edge_links) const {
  pos_edge_links->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Edge_range::const_iterator ei = mesh.edges().begin();
    ei != mesh.edges().end(); ++ei) {
    halfedge_descriptor hd = mesh.halfedge(*ei);
    if (mesh_properties->get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = mesh.opposite(hd);
    }
    const Link_list &edge_out_links = 
        mesh_properties->get_halfedge_out_links(hd);
    for (auto it = edge_out_links.begin();
      it != edge_out_links.end(); ++it) {
      const Link &link = *it;
      const Point &p = link.second.first;
      const Point &q = link.second.second;
      compute_segment(p, q, pos_edge_links);
    }
  }
}

void Scene::compute_vertex_start_points(bool is_input,
    std::vector<float> *pos_vertex_start_points) const {
  pos_vertex_start_points->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
    vi != mesh.vertices().end(); ++vi) {
    const Link &link = mesh_properties->get_vertex_out_link(*vi);
    compute_point(link.second.first, pos_vertex_start_points);
  }
}

void Scene::compute_vertex_end_points(bool is_input,
    std::vector<float> *pos_vertex_end_points) const {
  pos_vertex_end_points->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
    vi != mesh.vertices().end(); ++vi) {
    const Link &link = mesh_properties->get_vertex_out_link(*vi);
    compute_point(link.second.second, pos_vertex_end_points);
  }
}

void Scene::compute_vertex_links(bool is_input,
    std::vector<float> *pos_vertex_links) const {
  pos_vertex_links->clear();
  const Mesh_properties *mesh_properties = get_mesh_properties(is_input);
  const Mesh &mesh = mesh_properties->get_mesh();
  for (Mesh::Vertex_range::const_iterator vi = mesh.vertices().begin();
      vi != mesh.vertices().end(); ++vi) {
    const Link &link = mesh_properties->get_vertex_out_link(*vi);
    const Point &p = link.second.first;
    const Point &q = link.second.second;
    compute_segment(p, q, pos_vertex_links);
  }
}

void inline Scene::compute_segment(const Point &p, const Point &q,
  std::vector<float> *pos) const {
  compute_point(p, pos);
  compute_point(q, pos);
}

void inline Scene::compute_point(const Point &point,
  std::vector<float> *pos) const {
  pos->push_back(point.x());
  pos->push_back(point.y());
  pos->push_back(point.z());
}

void inline Scene::compute_normal(const Normal &normal,
  std::vector<float> *pos) const {
  pos->push_back(normal.x());
  pos->push_back(normal.y());
  pos->push_back(normal.z());
}

void inline Scene::compute_color(const Color &color,
  std::vector<float> *pos) const {
  pos->push_back(color.red() / 255.0f);
  pos->push_back(color.green() / 255.0f);
  pos->push_back(color.blue() / 255.0f);
}

Color Scene::get_rgb_color(FT h, FT s, FT v) const {
  // h(1~360), s(0~1), v(0~1)
  int hi = h / 60;
  FT f = h / 60 - hi, p = v * (1.0 - s);
  FT q = v * (1.0 - f * s), t = v * (1.0 - (1 - f) * s);
  v *= 255, t *= 255, p *= 255, q *= 255;
  switch (hi) {
  case 0:
    return Color(v, t, p);
  case 1:
    return Color(q, v, p);
  case 2:
    return Color(p, v, t);
  case 3:
    return Color(p, q, v);
  case 4:
    return Color(t, p, v);
  case 5:
    return Color(v, p, q);
  default:
    return Color(0, 0, 0);
  }
}

bool Scene::read_ply(std::ifstream &in, Mesh *mesh) const {
  /* Now we only support:
     the ply file based on ASCII code
     the triangles (no quad), and vertex_indices should be the 1st property
     all properties for one element is in one line
  */
  // precondition: in is valid

  std::set<halfedge_descriptor> crease_halfedges;
  std::string line;
  std::vector<std::string> tokens;

  // step 1: read the head
  std::getline(in, line);     // "ply"
  std::getline(in, line);     // "format ascii 1.0"
  boost::split(tokens, line, boost::is_any_of(" "));
  if (tokens[1] != "ascii") {
    return false;             // not in ASCII format
  }
  // type_name, type_count. e.g., vertex -> 8, face -> 7, edge -> 5
  std::vector<std::pair<std::string, int>> element_types;
  // count, type, name. for a scaler, count == 1; for a list, count = 4
  std::vector<std::vector<std::vector<std::string>>> element_properties;
  while (!in.eof()) {
    std::getline(in, line);
    boost::trim(line);
    if (line.length() == 0) {
      continue;
    }
    boost::split(tokens, line, boost::is_any_of(" "));
    if (tokens[0] == "comment") {
      continue;
    } else if (tokens[0] == "element") {
      // try to extract vertex, face and edge information
      element_types.push_back(std::make_pair(tokens[1], std::stoi(tokens[2])));
      element_properties.push_back({});
    } else if (tokens[0] == "end_header") {
      break;
    } else if (tokens[0] == "property") {
      std::vector<std::string> element_property;
      if (tokens[1] == "list") {  // property is a list
        element_property.push_back("4");
        element_property.push_back(tokens[3]);    // property data type
        element_property.push_back(tokens[4]);    // property name
      }
      else {                      // property is a scaler
        element_property.push_back("1");          // property count
        element_property.push_back(tokens[1]);    // property data type
        element_property.push_back(tokens[2]);    // property name
      }
      element_properties.back().push_back(element_property);
    }
  }

  // step 2: read the data
  std::vector<vertex_descriptor> vertex_descriptors;
  for (int i = 0; i < element_types.size(); ++i) {
    if (element_types[i].first == "vertex") {
      // extract the index of x, y and z
      int x_index = -1, y_index = -1, z_index = -1;
      int index = 0;
      for (int j = 0; j < element_properties[i].size(); ++j) {
        if (element_properties[i][j][2] == "x") {
          x_index = index;
        }
        else if (element_properties[i][j][2] == "y") {
          y_index = index;
        }
        else if (element_properties[i][j][2] == "z") {
          z_index = index;
        }
        index += std::stoi(element_properties[i][j][0]);
      }
      // read the data, construct and add the vertex
      for (int j = 0; j < element_types[i].second; ++j) {
        std::getline(in, line);
        boost::trim(line);
        boost::split(tokens, line, boost::is_any_of(" "));
        Point p(std::stod(tokens[x_index]),
                std::stod(tokens[y_index]),
                std::stod(tokens[z_index]));
        vertex_descriptor vd = mesh->add_vertex(p);
        vertex_descriptors.push_back(vd);
      }
    }
    else if (element_types[i].first == "face") {
      // extract the index of vertex_index
      int vertex_index = -1;
      int index = 0;
      for (int j = 0; j < element_properties[i].size(); ++j) {
        if (element_properties[i][j][2] == "vertex_index" ||
            element_properties[i][j][2] == "vertex_indices") {
          vertex_index = index;
          break;
        }
        index += std::stoi(element_properties[i][j][0]);
      }
      // read the data, construct and add the face
      for (int j = 0; j < element_types[i].second; ++j) {
        std::getline(in, line);
        boost::trim(line);
        boost::split(tokens, line, boost::is_any_of(" "));
        int vertice_count = std::stoi(tokens[vertex_index]);
        if (vertice_count != 3) { // we only prcess triangles here
          return false;
        }
        int index1 = std::stoi(tokens[vertex_index + 1]);
        int index2 = std::stoi(tokens[vertex_index + 2]);
        int index3 = std::stoi(tokens[vertex_index + 3]);
        mesh->add_face(vertex_descriptors[index1],
                       vertex_descriptors[index2], 
                       vertex_descriptors[index3]);
      }
    }
    else if (element_types[i].first == "edge") {
      // extract the index of vertex1, vertex2 and is_crease
      int vertex1_index = -1, vertex2_index = -1, is_crease_index = -1;
      int index = 0;
      for (int j = 0; j < element_properties[i].size(); ++j) {
        if (element_properties[i][j][2] == "vertex_1") {
          vertex1_index = index;
        }
        else if (element_properties[i][j][2] == "vertex_2") {
          vertex2_index = index;
        }
        else if (element_properties[i][j][2] == "is_crease") {
          is_crease_index = index;
        }
        index += std::stoi(element_properties[i][j][0]);
      }
      // read the data, add the crease edges to crease_halfedges
      for (int j = 0; j < element_types[i].second; ++j) {
        std::getline(in, line);
        boost::trim(line);
        boost::split(tokens, line, boost::is_any_of(" "));
        int v1_index = std::stoi(tokens[vertex1_index]);
        int v2_index = std::stoi(tokens[vertex2_index]);
        int is_crease = std::stoi(tokens[is_crease_index]);
        if (is_crease != 0) {
          vertex_descriptor vd1 = vertex_descriptors[v1_index];
          vertex_descriptor vd2 = vertex_descriptors[v2_index];
          Halfedge_around_target_circulator hb(mesh->halfedge(vd2), *mesh);
          Halfedge_around_target_circulator he(hb);
          do {
            if (mesh->source(*hb) == vd1) {
              crease_halfedges.insert(*hb);
              break;
            }
            ++hb;
          } while (hb != he);
        }
      }
    }
  }

  // step 3: add the feature halfedge property if necessary
  if (!crease_halfedges.empty()) {
    Mesh::Property_map<halfedge_descriptor, bool> halfedge_are_creases;
    bool created;
    boost::tie(halfedge_are_creases, created) =
        mesh->add_property_map<halfedge_descriptor, bool>("h:crease", false);
    assert(created);
    for (Mesh::Edge_range::const_iterator ei = mesh->edges().begin();
        ei != mesh->edges().end(); ++ei) {
      halfedge_descriptor hd = mesh->halfedge(*ei);
      if (!mesh->is_border(hd)) {
        hd = mesh->opposite(hd);
      }
      if (mesh->is_border(hd)) {  // ei is a boundary edge (add automatically)
        // add the boundary halfedge as crease automatically
        halfedge_are_creases[hd] = false;
        halfedge_are_creases[mesh->opposite(hd)] = true;
      } else {                    // ei is an inner edge (check in the set)
        // add the specified halfege as crease if found in crease_halfedges
        auto it = crease_halfedges.find(hd);
        if (it == crease_halfedges.end()) { // its opposite is recorded
          hd = mesh->opposite(hd);
          it = crease_halfedges.find(hd);
        }
        if (it != crease_halfedges.end()) {
          halfedge_are_creases[hd] = true;
          halfedge_are_creases[mesh->opposite(hd)] = false;
        }
      }
    }
  }
  return true;
  // TODO: understand what the function "add_property_map" does.
  //mesh->add_property_map();
  // please refer to https://doc.cgal.org/latest/Surface_mesh/index.html
}

bool Scene::open_surface_mesh(QString file_name, Mesh *mesh) const {
  QTextStream cerr(stderr);
  cerr << QString("\nOpening file \"%1\"\n").arg(file_name);
  cerr.flush();

  QFileInfo file_info(file_name);
  std::ifstream in(file_name.toUtf8());
  if (!file_info.isFile() || !file_info.isReadable() || !in) {
    std::cerr << "unable to open file" << std::endl;
    return false;
  }

  if (file_name.endsWith(".off", Qt::CaseInsensitive)) {
    in >> *mesh;
    if (!in) {
      std::cerr << "invalid OFF file" << std::endl;
      return false;
    }
    in.close();
  }
  else if (file_name.endsWith(".ply", Qt::CaseInsensitive)) {
    // read the ply file, and store the crease edges in property map "h:crease"
    bool suc = read_ply(in, mesh);
    in.close();
    if (!suc) {
      std::cerr << "invalid PLY file" << std::endl;
      return false;
    }
  }
  else {
    std::cerr << "file format not supported" << std::endl;
    in.close();
    return false;
  }
  return true;
}

bool Scene::open(QString file_name) {
  // step 1: open the file, construct the mesh and crease halfedges
  Mesh *mesh = new Mesh();
  bool suc = open_surface_mesh(file_name, mesh);
  if (!suc) {
    delete mesh;
    return false;
  }
  // step 2: set the m_pInput and m_pRemesh
  if (m_pInput != NULL) {
    m_minangle_remesh.delete_input();
    delete m_pInput;
  }
  if (m_pRemesh != NULL) {
    m_minangle_remesh.delete_remesh();
    delete m_pRemesh;
  }
  m_pInput = mesh;
  normalize(1.0, m_pInput);
  m_minangle_remesh.set_input(m_pInput, true);
  m_target_edge_length = calculate_input_edge_length();
  update_bbox();
  m_pRemesh = new Mesh(*m_pInput);
  m_minangle_remesh.set_remesh(m_pRemesh, false);
  reset_draw_render_types();
  m_view_input = false;
  m_view_remesh = true;

  changed();
  return true;
}

bool Scene::open_input(QString file_name) {
  // step 1: open the file, construct the mesh and crease halfedges
  Mesh *mesh = new Mesh();
  bool suc = open_surface_mesh(file_name, mesh);
  if (!suc) {
    delete mesh;
    return false;
  }
  // step 2: set the m_pInput, and m_pRemesh if necessary
  if (m_pInput != NULL) {
    m_minangle_remesh.delete_input();
    delete m_pInput;
  }
  m_pInput = mesh;
  normalize(1.0, m_pInput);
  m_minangle_remesh.set_input(m_pInput, true);
  m_target_edge_length = calculate_input_edge_length();
  update_bbox();
  if (m_pRemesh == NULL) {
    m_pRemesh = new Mesh(*m_pInput);
    m_minangle_remesh.set_remesh(m_pRemesh, false);
  }
  reset_draw_render_types();
  m_view_input = true;
  m_view_remesh = false;

  changed();
  return true;
}

bool Scene::open_remesh(QString file_name) {
  // step 1: open the file, construct the mesh and crease_halfedges
  Mesh *mesh = new Mesh();
  bool suc = open_surface_mesh(file_name, mesh);
  if (!suc) {
    delete mesh;
    return false;
  }
  // step 2: set the m_pRemesh, and m_pInput if necessary
  if (m_pRemesh != NULL) {
    m_minangle_remesh.delete_remesh();
    delete m_pRemesh;
  }
  m_pRemesh = mesh;
  if (m_pInput == NULL) {
    normalize(1.0, m_pRemesh);
    m_pInput = new Mesh(*m_pRemesh);
    m_minangle_remesh.set_input(m_pInput, false);
    m_target_edge_length = calculate_input_edge_length();
    update_bbox();
  }
  m_minangle_remesh.set_remesh(m_pRemesh, true);
  reset_draw_render_types();
  m_view_input = false;
  m_view_remesh = true;

  changed();
  return true;
}

void Scene::save_remesh_as(QString file_name) {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.save_remesh_as(file_name.toStdString());    
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

void Scene::toggle_view_mesh_edges() {
  m_view_mesh_edges = !m_view_mesh_edges;
  if (m_view_mesh_edges) {
    changed();
  }
}

void Scene::toggle_view_mesh_plain_faces() {
  reset_draw_render_types();
  changed();
}

void Scene::toggle_view_interpolated_feature_intensities() {
  set_draw_render_types(DrawType::k_mesh, RenderType::k_ifi_faces);
  changed();
}

void Scene::toggle_view_face_errors() {
  set_draw_render_types(DrawType::k_mesh, RenderType::k_mr_faces);
  changed();
}

void Scene::toggle_view_element_classifications() {
  set_draw_render_types(DrawType::k_vertex_voronoi, 
                        RenderType::k_classifications);
  changed();
}

void Scene::toggle_view_gaussian_curvatures() {
  set_draw_render_types(DrawType::k_vertex_voronoi,
                        RenderType::k_gaussian_curvature);
  changed();
}

void Scene::toggle_view_maximal_normal_dihedrals() {
  set_draw_render_types(DrawType::k_vertex_voronoi,
                        RenderType::k_maximal_halfedge_dihedral);
  changed();
}

void Scene::toggle_view_normal_dihedrals() {
  set_draw_render_types(DrawType::k_edge_voronoi,
                        RenderType::k_normal_dihedral);
  changed();
}

void Scene::toggle_view_face_in_start_points() {
  m_view_face_in_start_points = !m_view_face_in_start_points;
  if (m_view_face_in_start_points) {
    changed();
  }
}

void Scene::toggle_view_face_in_end_points() {
  m_view_face_in_end_points = !m_view_face_in_end_points;
  if (m_view_face_in_end_points) {
    changed();
  }
}

void Scene::toggle_view_face_in_links() {
  m_view_face_in_links = !m_view_face_in_links;
  if (m_view_face_in_links) {
    changed();
  }
  changed();
}

void Scene::toggle_view_face_out_start_points() {
  m_view_face_out_start_points = !m_view_face_out_start_points;
  if (m_view_face_out_start_points) {
    changed();
  }
}

void Scene::toggle_view_face_out_end_points() {
  m_view_face_out_end_points = !m_view_face_out_end_points;
  if (m_view_face_out_end_points) {
    changed();
  }
}

void Scene::toggle_view_face_out_links() {
  m_view_face_out_links = !m_view_face_out_links;
  if (m_view_face_out_links) {
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

void Scene::toggle_view_face_feature_intensities() {      // face samples
  set_draw_render_types(k_face_voronoi, k_feature_intensity);
  changed();
}

void Scene::toggle_view_face_capacities() {
  set_draw_render_types(k_face_voronoi, k_capacity);
  changed();
}

void Scene::toggle_view_face_weights() {
  set_draw_render_types(k_face_voronoi, k_weight);
  changed();
}

void Scene::eliminate_degenerations() {
  if (m_pInput == NULL) {
    std::cout << "Please open a file first" << std::endl;
    return;
  }
  CGAL::Timer timer;
  timer.start();
  std::cout << std::endl << "Eliminate Input degenerated faces...";
  Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
  size_t nb_eliminations = remesher->eliminate_input_degenerated_faces();
  std::cout << "Done (" << nb_eliminations << " faces eliminated, "
    << timer.time() << " s)" << std::endl;
  if (nb_eliminations > 0) {
    m_minangle_remesh.set_input(m_pInput, false);
    m_target_edge_length = calculate_input_edge_length();
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
  Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
  size_t nb_split = remesher->split_input_long_edges();
  std::cout << "Done (" << nb_split << " edges splited, "
    << timer.time() << " s)" << std::endl;
  if (nb_split > 0) {
    m_minangle_remesh.set_input(m_pInput, false);
    m_target_edge_length = calculate_input_edge_length();
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
    m_minangle_remesh.get_remesher()->input_properties();
  }
}

void Scene::split_borders() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Splitting borders...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(
        faces(*m_pRemesh),
        *m_pRemesh,
        boost::make_function_output_iterator(halfedge2edge(*m_pRemesh, border)));
    PMP::split_long_edges(border, m_target_edge_length, *m_pRemesh);
    m_minangle_remesh.set_remesh(m_pRemesh, false);
    std::cout << "Done (" << timer.time() << " s)" << std::endl;
    reset_draw_render_types();
    changed();
  }
}

void Scene::isotropic_remeshing() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Isotropic remeshing...";
    PMP::isotropic_remeshing(
        faces(*m_pRemesh), 
        m_target_edge_length, 
        *m_pRemesh,
        PMP::parameters::number_of_iterations(m_smooth_iteration_count)
        .protect_constraints(true)  //i.e. protect border, here
      );
    m_minangle_remesh.set_remesh(m_pRemesh, false);
    std::cout << "Done (" << timer.time() << " s)" << std::endl;
    reset_draw_render_types();
    changed();
  }
}

void Scene::test() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Minimal angle remeshing...";
    /* PMP::minimal_angle_remeshing(*m_pRemesh, 
      CGAL::parameters::min_angle_threshold(15.0).verbose_progress(false)); */
    PMP::minimal_angle_remeshing(*m_pRemesh);
    m_minangle_remesh.set_remesh(m_pRemesh, false);
    std::cout << "Done (" << timer.time() << " s)" << std::endl;
    reset_draw_render_types();
    changed();
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

void Scene::generate_links() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
    return;
  }
  else {
    m_minangle_remesh.generate_samples_and_links();
    changed();
  }
}

void Scene::remesh_properties() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    bool new_generated = m_minangle_remesh.get_remesher()->remesh_properties();
    if (new_generated) {
      changed();
    }
  }
}

void Scene::minangle_remeshing() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.minangle_remeshing();
    reset_draw_render_types();
    changed();
  }
}

void Scene::initial_mesh_simplification() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.initial_mesh_simplification();
    reset_draw_render_types();
    changed();
  }
}

void Scene::split_local_longest_edge() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.get_remesher()->split_local_longest_edge();
    reset_draw_render_types();
    changed();
  }
}

void Scene::increase_minimal_angle() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.increase_minimal_angle();
    reset_draw_render_types();
    changed();
  }
}

void Scene::maximize_minimal_angle() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.maximize_minimal_angle();
    reset_draw_render_types();
    changed();
  }
}

void Scene::final_vertex_relocation() {
  if (m_pRemesh == NULL) {
    std::cout << "Please open a file first" << std::endl;
  }
  else {
    m_minangle_remesh.final_vertex_relocation();
    reset_draw_render_types();
    changed();
  }
}

void Scene::update_feature_intensities() {
  // links will be cleared automatically once feature intensities updated
  if (m_pInput != NULL && m_pRemesh != NULL) {
    Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
    remesher->calculate_feature_intensities(true, true, true);
    changed();
  }
}

void Scene::set_draw_render_types(DrawType draw_type, RenderType render_type) {
  m_draw_type = draw_type;
  m_render_type = render_type;
}

void Scene::reset_draw_render_types() {
  set_draw_render_types(k_mesh, k_plain_faces);
}

void Scene::reset() {
  if (m_pInput != NULL) {
    if (m_pRemesh != NULL) {
      delete m_pRemesh;
    }
    m_pRemesh = new Mesh(*m_pInput);
    m_minangle_remesh.set_remesh(m_pRemesh, false);
  }
}

void Scene::normalize(FT radius, Mesh *mesh) const {
  // step 1: calculate the bounding box
  Bbox bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
                   DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
  Mesh::Vertex_range::const_iterator vi = mesh->vertices().begin();
  bbox = mesh->point(*vi).bbox();
  for (; vi != mesh->vertices().end(); ++vi) {
    bbox = bbox + mesh->point(*vi).bbox();
  }
  // step 2: get the center and radius
  FT x_center = (bbox.xmin() + bbox.xmax()) / 2.0;
  FT y_center = (bbox.ymin() + bbox.ymax()) / 2.0;
  FT z_center = (bbox.zmin() + bbox.zmax()) / 2.0;
  FT x_radius = (bbox.xmax() - bbox.xmin()) / 2.0;
  FT y_radius = (bbox.ymax() - bbox.ymin()) / 2.0;
  FT z_radius = (bbox.zmax() - bbox.zmin()) / 2.0;
  FT max_radius = std::max(std::max(x_radius, y_radius), z_radius);
  // step 3: transfer
  for (Mesh::Vertex_range::const_iterator vi = mesh->vertices().begin();
    vi != mesh->vertices().end(); ++vi) {
    Point &p = mesh->point(*vi);
    p = Point(p.x() - x_center, p.y() - y_center, p.z() - z_center);
  }
  // step 4: scale
  for (Mesh::Vertex_range::const_iterator vi = mesh->vertices().begin();
    vi != mesh->vertices().end(); ++vi) {
    Point &p = mesh->point(*vi);
    p = Point(p.x() * radius / max_radius,
              p.y() * radius / max_radius,
              p.z() * radius / max_radius);
  }
}

double Scene::calculate_input_edge_length() const {
  const Mesh_properties *mesh_properties = get_mesh_properties(true);
  if (mesh_properties == NULL) {
    return TARGET_EDGE_LENGTH;
  }
  else {
    return mesh_properties->calculate_average_length();
  }
}

