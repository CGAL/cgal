#include "Viewer.h"

#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>

#include <CGAL/Qt/CreateOpenGLContext.h>
#include "create_sphere.h"

#include <vector>

Viewer::Viewer(QWidget* parent )
  : CGAL::QGLViewer(parent)
{
  draw_balls = true;
  draw_inputs = false;
  trivial_center.push_back(0.);
  trivial_center.push_back(0.);
  trivial_center.push_back(0.);
  radius_ = 1.0;
  center_ = Point_3(0,0,0);
}

void Viewer::compile_shaders()
{
  initializeOpenGLFunctions();
  for(int i = 0; i<SIZE_OF_VBO; ++i)
    if(! buffers[i].create())
    {
      std::cerr<<"VBO Creation FAILED"<<std::endl;
      return;
    }
  for(int i = 0; i<SIZE_OF_VAO; ++i)
    if(!vao[i].create())
    {
      std::cerr<<"VAO Creation FAILED"<<std::endl;
      return;
    }

  //The sphere

  //Vertex source code
  const char vertex_source[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 normal;\n"
    "attribute highp vec4 center;\n"

    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix; \n"

    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "void main(void)\n"
    "{\n"
    "   fP = mv_matrix * (vertex+vec4(center.xyz, 0.0)); \n"
    "   fN = mat3(mv_matrix)* normal; \n"
    "   gl_Position = mvp_matrix * (vertex+vec4(center.xyz, 0.0));\n"
    "}"
  };
  //Vertex source code
  const char fragment_source[] =
  {
    "#version 120 \n"
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "uniform highp vec4 color; \n"
    "uniform highp vec4 light_pos;  \n"
    "uniform highp vec4 light_diff; \n"
    "uniform highp vec4 light_spec; \n"
    "uniform highp vec4 light_amb;  \n"
    "uniform float spec_power ; \n"

    "void main(void) { \n"

    "   highp vec3 L = light_pos.xyz - fP.xyz; \n"
    "   highp vec3 V = -fP.xyz; \n"

    "   highp vec3 N = normalize(fN); \n"
    "   L = normalize(L); \n"
    "   V = normalize(V); \n"

    "   highp vec3 R = reflect(-L, N); \n"
    "   highp vec4 diffuse = max(dot(N,L), 0.0) * light_diff * color; \n"
    "   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

    "gl_FragColor = light_amb*color + diffuse  ; \n"
    "} \n"
    "\n"
  };
  QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
  if(!vertex_shader->compileSourceCode(vertex_source))
  {
    std::cerr<<"Compiling vertex source FAILED"<<std::endl;
  }

  QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
  if(!fragment_shader->compileSourceCode(fragment_source))
  {
    std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
  }

  if(!rendering_program.addShader(vertex_shader))
  {
    std::cerr<<"adding vertex shader FAILED"<<std::endl;
  }
  if(!rendering_program.addShader(fragment_shader))
  {
    std::cerr<<"adding fragment shader FAILED"<<std::endl;
  }
  if(!rendering_program.link())
  {
    std::cerr<<"linking Program FAILED"<<std::endl;
  }

      //The edges

      //Vertex source code
      const char vertex_source_edges[] =
      {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * vertex;\n"
        "}"
      };
      //Vertex source code
      const char fragment_source_edges[] =
      {
        "#version 120 \n"
        "void main(void) { \n"
        "gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0); \n"
        "} \n"
        "\n"
      };
      QOpenGLShader *vertex_shader_edges = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_edges->compileSourceCode(vertex_source_edges))
      {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }

      QOpenGLShader *fragment_shader_edges = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_edges->compileSourceCode(fragment_source_edges))
      {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }

      if(!rendering_program_edges.addShader(vertex_shader_edges))
      {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
      }
      if(!rendering_program_edges.addShader(fragment_shader_edges))
      {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
      }
      if(!rendering_program_edges.link())
      {
        std::cerr<<"linking Program FAILED"<<std::endl;
      }
}

void Viewer::initialize_buffers()
{
  rendering_program.bind();
  //The big white sphere
  vao[SPHERE_VAO].bind();

  //points of the sphere
  buffers[SPHERE_POINTS].bind();
  buffers[SPHERE_POINTS].allocate(pos_sphere.data(),
                                  static_cast<int>(pos_sphere.size()*sizeof(float)));
  vertexLocation[0] = rendering_program.attributeLocation("vertex");
  rendering_program.enableAttributeArray(vertexLocation[0]);
  rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
  buffers[SPHERE_POINTS].release();

  //normals of the sphere
  buffers[SPHERE_NORMALS].bind();
  buffers[SPHERE_NORMALS].allocate(normals.data(),
                                   static_cast<int>(normals.size()*sizeof(float)));
  normalsLocation[0] = rendering_program.attributeLocation("normal");
  rendering_program.enableAttributeArray(normalsLocation[0]);
  rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
  buffers[SPHERE_NORMALS].release();

  //center of the sphere
  buffers[SPHERE_CENTER].bind();
  buffers[SPHERE_CENTER].allocate(trivial_center.data(),
                                  static_cast<int>(trivial_center.size()*sizeof(float)));
  trivialCenterLocation = rendering_program.attributeLocation("center");
  rendering_program.enableAttributeArray(trivialCenterLocation);
  rendering_program.setAttributeBuffer(trivialCenterLocation,GL_FLOAT,0,3);
  buffers[SPHERE_CENTER].release();

  glVertexAttribDivisor(trivialCenterLocation, 1);
  glVertexAttribDivisor(normalsLocation[0], 0);
  vao[SPHERE_VAO].release();

  //The little green spheres
  vao[POINTS_VAO].bind();

  //points of the spheres
  buffers[POINTS_POINTS].bind();
  buffers[POINTS_POINTS].allocate(pos_sphere_inter.data(),
                                  static_cast<int>(pos_sphere_inter.size()*sizeof(float)));
  vertexLocation[2] = rendering_program.attributeLocation("vertex");
  rendering_program.enableAttributeArray(vertexLocation[2]);
  rendering_program.setAttributeBuffer(vertexLocation[2],GL_FLOAT,0,3);
  buffers[POINTS_POINTS].release();

  //normals of the sphere
  buffers[POINTS_NORMALS].bind();
  buffers[POINTS_NORMALS].allocate(normals_inter.data(),
                                   static_cast<int>(normals_inter.size()*sizeof(float)));
  normalsLocation[2] = rendering_program.attributeLocation("normal");
  rendering_program.enableAttributeArray(normalsLocation[2]);
  rendering_program.setAttributeBuffer(normalsLocation[2],GL_FLOAT,0,3);
  buffers[POINTS_NORMALS].release();

  //center of the sphere
  buffers[POINTS_CENTERS].bind();
  buffers[POINTS_CENTERS].allocate(pos_points.data(),
                                   static_cast<int>(pos_points.size()*sizeof(float)));
  centerLocation = rendering_program.attributeLocation("center");
  rendering_program.enableAttributeArray(centerLocation);
  rendering_program.setAttributeBuffer(centerLocation,GL_FLOAT,0,3);
  buffers[POINTS_CENTERS].release();

  glVertexAttribDivisor(centerLocation, 1);
  glVertexAttribDivisor(normalsLocation[2], 0);
  vao[POINTS_VAO].release();
  rendering_program.release();

  //The circles
  rendering_program_edges.bind();
  vao[EDGES_VAO].bind();
  buffers[EDGES_POINTS].bind();
  buffers[EDGES_POINTS].allocate(pos_lines.data(),
                                 static_cast<int>(pos_lines.size()*sizeof(float)));
  vertexLocation[1] = rendering_program_edges.attributeLocation("vertex");
  rendering_program_edges.enableAttributeArray(vertexLocation[1]);
  rendering_program_edges.setAttributeBuffer(vertexLocation[1],GL_FLOAT,0,3);
  buffers[EDGES_POINTS].release();
  vao[EDGES_VAO].release();
  rendering_program_edges.release();
}

void Viewer::compute_elements()
{
  pos_sphere_inter.clear();
  normals_inter.clear();
  pos_sphere.clear();
  normals.clear();
  trivial_center.clear();
  trivial_center.push_back(static_cast<float>(center_.x()));
  trivial_center.push_back(static_cast<float>(center_.y()));
  trivial_center.push_back(static_cast<float>(center_.z()));
  create_flat_sphere(static_cast<float>(radius_), pos_sphere, normals, 9);
  create_flat_sphere(0.001f*static_cast<float>(radius_), pos_sphere_inter, normals_inter, 36);

  pos_points.resize(0);
  pos_lines.resize(0);

  std::cout << inputs.size() << " inputs to points" << std::endl;
  // draw points as small spheres
  for(const auto& p : inputs)
  {
    pos_points.push_back(p.x());
    pos_points.push_back(p.y());
    pos_points.push_back(p.z());
  }

  std::cout << subsampled_arcs.size() << " edge arcs" << std::endl;

  // draw the primal edges
  std::size_t tot_size = 0;
  for(const auto& arc : subsampled_arcs)
    tot_size += arc.size();
  pos_lines.reserve(2 * 3 * tot_size);

  for(const auto& arc : subsampled_arcs)
  {
    for(std::size_t i=0; i<arc.size()-1; ++i)
    {
      pos_lines.push_back(arc[i].x());
      pos_lines.push_back(arc[i].y());
      pos_lines.push_back(arc[i].z());
      pos_lines.push_back(arc[i+1].x());
      pos_lines.push_back(arc[i+1].y());
      pos_lines.push_back(arc[i+1].z());
    }
  }
  setSceneRadius(radius_*1.3);
  setSceneCenter({center_.x(), center_.y(), center_.z()});
  showEntireScene();
}

void Viewer::attrib_buffers(CGAL::QGLViewer* viewer)
{
  QMatrix4x4 mvpMatrix;
  QMatrix4x4 mvMatrix;
  double mat[16];

  viewer->camera()->getModelViewProjectionMatrix(mat);
  for(int i=0; i < 16; i++)
    mvpMatrix.data()[i] = (float)mat[i];

  viewer->camera()->getModelViewMatrix(mat);

  for(int i=0; i < 16; i++)
    mvMatrix.data()[i] = (float)mat[i];

  // define material
  QVector4D        ambient(0.1f, 0.1f, 0.1f, 1.0f);
  QVector4D        diffuse( 0.9f,
                            0.9f,
                            0.9f,
                            0.0f );

  QVector4D        specular(  0.0f,
                              0.0f,
                              0.0f,
                              0.0f );

  QVector4D        position( -1.2f, 1.2f, .9797958971f, 1.0f  );
  GLfloat shininess =  1.0f;


  mvpLocation = rendering_program.uniformLocation("mvp_matrix");
  mvLocation = rendering_program.uniformLocation("mv_matrix");
  colorLocation = rendering_program.uniformLocation("color");
  lightLocation[0] = rendering_program.uniformLocation("light_pos");
  lightLocation[1] = rendering_program.uniformLocation("light_diff");
  lightLocation[2] = rendering_program.uniformLocation("light_spec");
  lightLocation[3] = rendering_program.uniformLocation("light_amb");
  lightLocation[4] = rendering_program.uniformLocation("spec_power");

  rendering_program.setUniformValue(lightLocation[0], position);
  rendering_program.setUniformValue(lightLocation[1], diffuse);
  rendering_program.setUniformValue(lightLocation[2], specular);
  rendering_program.setUniformValue(lightLocation[3], ambient);
  rendering_program.setUniformValue(lightLocation[4], shininess);
  rendering_program.setUniformValue(mvpLocation, mvpMatrix);
  rendering_program.setUniformValue(mvLocation, mvMatrix);
  rendering_program_edges.setUniformValue("mvp_matrix", mvpMatrix);

}

void Viewer::draw()
{
  glEnable(GL_DEPTH_TEST);
  QColor color;
  rendering_program.bind();

  if(draw_balls)
  {
    vao[SPHERE_VAO].bind();
    attrib_buffers(this);
    color.setRgbF(1.0f, 1.0f, 1.0f);
    rendering_program.setUniformValue(colorLocation, color);
    glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_sphere.size()/3), 1);
    vao[SPHERE_VAO].release();
  }

  if(draw_inputs)
  {
    vao[POINTS_VAO].bind();
    attrib_buffers(this);
    color.setRgbF(0.0f, 1.0f, 0.0f);
    rendering_program.setUniformValue(colorLocation, color);

    glDrawArraysInstanced(GL_TRIANGLES, 0,
                          static_cast<GLsizei>(pos_sphere_inter.size()/3),
                          static_cast<GLsizei>(pos_points.size()/3));
    vao[POINTS_VAO].release();
  }
  rendering_program.release();
  rendering_program_edges.bind();
  vao[EDGES_VAO].bind();
  attrib_buffers(this);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
  vao[EDGES_VAO].release();
  rendering_program_edges.release();
}

void Viewer::init()
{
  glDrawArraysInstanced = (PFNGLDRAWARRAYSINSTANCEDARBPROC)this->context()->getProcAddress("glDrawArraysInstancedARB");
  glVertexAttribDivisor = (PFNGLVERTEXATTRIBDIVISORARBPROC)this->context()->getProcAddress("glVertexAttribDivisorARB");
  if(!glDrawArraysInstanced || !glDrawArraysInstanced)
  {
    qDebug()<<"Missing OpenGL extensions. Demo not available.";
    return;
  }

  compile_shaders();
  compute_elements();

  initialize_buffers();
  glEnable(GL_BLEND);

  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
}

QString Viewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
  if(e->key() == Qt::Key_P)
  {
    draw_inputs = !draw_inputs;
    update();
    return;
  }

  if(e->key() == Qt::Key_S)
  {
    draw_balls = !draw_balls;
    update();
    return;
  }

  QGLViewer::keyPressEvent(e);
}
