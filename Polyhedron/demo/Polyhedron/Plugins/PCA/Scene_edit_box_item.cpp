#include "Scene_edit_box_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLShaderProgram>
#include <QCursor>

using namespace CGAL::Three;
struct Scene_edit_box_item::vertex{
  int id;
  double *x;
  double *y;
  double *z;

  CGAL::Point_3<Kernel> position()const
  {
    return  CGAL::Point_3<Kernel>(*x,*y,*z);
  }
  double operator[](int i)
  {
    switch(i)
    {
    case 0:
      return *x;
    case 1:
      return *y;
    case 2:
      return *z;
    default:
      return 0;
    }
  }
};
struct Scene_edit_box_item::edge{
  vertex* source;
  vertex* target;
};
struct Scene_edit_box_item::face{
  vertex* vertices[4];
};

struct Scene_edit_box_item_priv{
  typedef CGAL::Simple_cartesian<double>  Kernel;
  enum VAOs{
    Edges = 0,
    Spheres,
    Faces,
    S_Edges,
    S_Spheres,
    S_Faces,
    P_Edges,
    P_Spheres,
    P_Faces,
    NumberOfVaos
  };
  enum VBOs{
    VertexEdges = 0,
    ColorsEdges,
    VertexSpheres,
    NormalSpheres,
    CenterSpheres,
    ColorsSpheres,
    VertexFaces,
    NormalFaces,
    ColorsFaces,
    S_Vertex,
    S_Normal,
    NumberOfVbos
  };
  enum HL_Primitive{
    VERTEX=0,
    EDGE,
    FACE,
    NO_TYPE
  };
  Scene_edit_box_item_priv(const Scene_interface *scene_interface, Scene_edit_box_item* ebi)
  {
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    ready_to_hl = true;
    scene = scene_interface;
    item = ebi;
    selection_on = false;
    Scene_item::Bbox bb = scene->bbox();
    double x=(bb.xmin()+bb.xmax())/2;
    double y=(bb.ymin()+bb.ymax())/2;
    double z=(bb.zmin()+bb.zmax())/2;
    center_ = qglviewer::Vec(x,y,z);
    relative_center_ = qglviewer::Vec(0,0,0);
    remodel_frame = new Scene_item::ManipulatedFrame();
    remodel_frame->setTranslationSensitivity(1.0);
    frame = new Scene_item::ManipulatedFrame();
    frame->setPosition(center_+offset);
    frame->setSpinningSensitivity(100.0); //forbid spinning
    constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
    constraint.setRotationConstraintDirection(qglviewer::Vec(.0,.0,.1));
    frame->setConstraint(&constraint);
    //create the sphere model
    pool[0] = bb.xmin(); pool[3] = bb.xmax();
    pool[1] = bb.ymin(); pool[4] = bb.ymax();
    pool[2] = bb.zmin(); pool[5] = bb.zmax();

    vertex_spheres.resize(0);
    normal_spheres.resize(0);
    create_flat_sphere(1.0f, vertex_spheres, normal_spheres,10);
    //      5-----6
    //  .   |  .  |
    // 4------7   |
    // |    | |   |
    // |    1-|---2
    // | .    |.
    // 0------3

    //vertices
    for( int i = 0; i< 8; ++i)
    {

      vertices[i].x = ((i/2)%2==0)? &pool[0]:&pool[3];
      vertices[i].y = (i/4==0)? &pool[1]:&pool[4];
      vertices[i].z = (((i+1)/2)%2==1)? &pool[2]:&pool[5];
      vertices[i].id = i;
    }

    //      .--5--.
    //  4   |  6  |
    // .--7-1-.   2
    // |    | |   |
    // 0    .-39--.
    // | 8    |10
    // .--11--.

    //edges
    for( int i=0; i<12; ++i)
    {
      if(i<4)
      {
        edges[i].source = &vertices[i];
        edges[i].target = &vertices[i+4];
      }
      else if(i<8)
      {
        edges[i].source = &vertices[i];
        edges[i].target = &vertices[(i+1)%4 +4];
      }
      else
      {
        edges[i].source = &vertices[i%4];
        edges[i].target = &vertices[(i+1) %4];
      }
      vertex_edges.resize(0);
    }

    //
    //                          5------6
    //                          |      |
    //      .-----.             |  2   |
    //  .   |5 .  |             |      |
    // .------.2  |      5------1------2------6------5
    // | 1  | |   |      |      |      |      |      |
    // |   4.-|-3-.      |   1  |   0  |  3   |   5  |
    // | .  0 |.         |      |      |      |      |
    // .------.          4------0------3------7------4
    //                          |      |
    //                          |  4   |
    //                          |      |
    //                          4------7


    //faces
    for( int i=0; i<4; ++i)
    {
      faces[0].vertices[i] = &vertices[i];
    }

    for( int i=1; i<4; ++i)
    {
      faces[i].vertices[0] = &vertices[i];
      faces[i].vertices[1] = &vertices[(i-1)];
      faces[i].vertices[2] = &vertices[i+3];
      faces[i].vertices[3] = &vertices[i+4];
    }

    faces[4].vertices[0] = &vertices[0];
    faces[4].vertices[1] = &vertices[3];
    faces[4].vertices[2] = &vertices[7];
    faces[4].vertices[3] = &vertices[4];

    for( int i=0; i<4; ++i)
    {
      faces[5].vertices[i] = &vertices[i+4];
    }

    vertex_faces.resize(0);
    normal_faces.resize(0);

    for(int i=0; i<8; ++i)
      for(int j=0; j<3; ++j)
        last_pool[i][j] = vertices[i][j];

    pick_sphere_program.addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_spheres.v");
    pick_sphere_program.addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_without_light.f");
    pick_sphere_program.bindAttributeLocation("colors", 1);
    pick_sphere_program.link();


    //Vertex source code
    const char vertex_source[] =
    {
      "#version 120 \n                         "
      "attribute highp vec4 vertex;            "
      "attribute highp vec3 normals;           "
      "attribute highp vec4 colors;            "
      "uniform highp mat4 mvp_matrix;          "
      "uniform highp mat4 mv_matrix;           "
      "varying highp vec4 fP;                  "
      "varying highp vec3 fN;                  "
      "varying highp vec4 color;               "
      "void main(void)                         "
      "{                                       "
      "   color = colors;                       "
      "   fP = mv_matrix * vertex;             "
      "   fN = mat3(mv_matrix)* normals;       "
      "   gl_Position = mvp_matrix * vertex;   "
      "}\n                                     "
      "\n                                      "
    };

    //Fragment source code
    const char fragment_source[] =
    {
      "#version 120 \n"
      "varying highp vec4 color;"
      "varying highp vec4 fP; "
      "varying highp vec3 fN; "
      "uniform highp vec4 light_pos;  "
      "uniform highp vec4 light_diff; "
      "uniform highp vec4 light_spec; "
      "uniform highp vec4 light_amb;  "
      "uniform highp float spec_power ; "
      "uniform int is_two_side; "
      "uniform bool is_selected;"
      "void main(void) {"
      "highp vec3 L = light_pos.xyz - fP.xyz;"
      "highp vec3 V = -fP.xyz;"
      "highp vec3 N;"
      "if(fN == highp vec3(0.0,0.0,0.0)) "
      "N = highp vec3(0.0,0.0,0.0);"
      "else "
      "N = normalize(fN);"
      "L = normalize(L);"
      "V = normalize(V);"
      "highp vec3 R = reflect(-L, N);"
      "vec4 diffuse;"
      "if(is_two_side == 1) "
      "diffuse = abs(dot(N,L)) * light_diff * color;"
      "else "
      "diffuse = max(dot(N,L), 0.0) * light_diff * color;"
      "highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;"
      "vec4 ret_color = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1);"
      "if(is_selected) "
      "gl_FragColor = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, color.a);"
      "else "
      "gl_FragColor = vec4(ret_color.rgb, color.a); }\n"
      "\n"
    };
    transparent_face_program.addShaderFromSourceCode(QOpenGLShader::Vertex,vertex_source);
    transparent_face_program.addShaderFromSourceCode(QOpenGLShader::Fragment,fragment_source);
    transparent_face_program.bindAttributeLocation("colors", 1);
    transparent_face_program.link();
    reset_selection();
    last_picked_id = -1;
    last_picked_type = -1;
    QPixmap pix(":/cgal/cursors/resources/rotate_around_cursor.png");
    rotate_cursor = QCursor(pix);
  }
  ~Scene_edit_box_item_priv(){
    delete frame;
    delete remodel_frame;
  }
  mutable std::vector<float> vertex_edges;
  mutable std::vector<float> color_edges;
  mutable std::vector<float> vertex_spheres;
  mutable std::vector<float> normal_spheres;
  mutable std::vector<float> center_spheres;
  mutable std::vector<float> color_spheres;
  mutable std::vector<float> vertex_faces;
  mutable std::vector<float> normal_faces;
  mutable std::vector<float> color_faces;
  mutable std::vector<float> hl_vertex;
  mutable std::vector<float> hl_normal;

  double pool[6];
  //id|coord
  double last_pool[8][3];
  bool ready_to_hl;


  qglviewer::ManipulatedFrame* frame;
  qglviewer::ManipulatedFrame* remodel_frame;
  qglviewer::Vec rf_last_pos;
  qglviewer::LocalConstraint constraint;
  qglviewer::Vec center_;
  qglviewer::Vec relative_center_;

  mutable QOpenGLShaderProgram pick_sphere_program;
  mutable QOpenGLShaderProgram transparent_face_program;
  mutable Scene_edit_box_item::vertex vertices[8];
  mutable Scene_edit_box_item::edge edges[12];
  mutable Scene_edit_box_item::face faces[6];

  std::vector<Scene_edit_box_item::vertex*> selected_vertices;
  void reset_selection();
  bool selection_on;
  void picking(int& type, int& id, Viewer_interface *viewer);

  mutable QOpenGLShaderProgram* program;
  void initializeBuffers(Viewer_interface *viewer)const;

  void computeElements() const;
  void draw_picking(Viewer_interface*);
  void remodel_box(const QVector3D &dir);
  double applyX(int id, double x, double dirx);
  double applyY(int id, double y, double diry);
  double applyZ(int id, double z, double dirz);
  const Scene_interface* scene;
  Scene_edit_box_item* item;
  QPoint picked_pixel;
  HL_Primitive hl_type;
  int last_picked_id;
  int last_picked_type;
  QCursor rotate_cursor;

};


Scene_edit_box_item::Scene_edit_box_item(const Scene_interface *scene_interface)
  :  Scene_item(Scene_edit_box_item_priv::NumberOfVbos,Scene_edit_box_item_priv::NumberOfVaos)

{
  d = new Scene_edit_box_item_priv(scene_interface, this);

  are_buffers_filled = false;
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->setMouseTracking(true);
}
QString Scene_edit_box_item::toolTip() const {

  return QString();
}

void Scene_edit_box_item::drawSpheres(Viewer_interface *viewer, const QMatrix4x4 f_matrix ) const
{
  vaos[Scene_edit_box_item_priv::Spheres]->bind();
  GLdouble d_mat[16];
  QMatrix4x4 mvp_mat;
  viewer->camera()->getModelViewProjectionMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mvp_mat.data()[i] = GLfloat(d_mat[i]);
  mvp_mat = mvp_mat*f_matrix;
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  QVector4D light_pos(0.0f,0.0f,1.0f, 1.0f );
  light_pos = light_pos*f_matrix;
  double radius =std::sqrt(
      (point(6,0) - point(0,0)) * (point(6,0) - point(0,0)) +
      (point(6,1) - point(0,1)) * (point(6,1) - point(0,1)) +
      (point(6,2) - point(0,2)) * (point(6,2) - point(0,2))) *0.02 ;

  attribBuffers(viewer, PROGRAM_SPHERES);
  d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
  d->program->bind();
  d->program->setUniformValue("mvp_matrix", mvp_mat);
  d->program->setUniformValue("mv_matrix", mv_mat);
  d->program->setUniformValue("light_pos", light_pos);
  d->program->setAttributeValue("radius",radius);
  d->program->setAttributeValue("colors", QColor(Qt::red));
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(d->vertex_spheres.size()/3),
                                static_cast<GLsizei>(8));
  d->program->release();
  vaos[Scene_edit_box_item_priv::Spheres]->release();
}

void Scene_edit_box_item::draw(Viewer_interface *viewer) const
{
  if (!are_buffers_filled)
  {
    d->computeElements();
    d->initializeBuffers(viewer);
  }
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }
  GLdouble d_mat[16];
  QMatrix4x4 mvp_mat;
  viewer->camera()->getModelViewProjectionMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mvp_mat.data()[i] = GLfloat(d_mat[i]);
  mvp_mat = mvp_mat*f_matrix;
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  QVector4D light_pos(0.0f,0.0f,1.0f, 1.0f );
  light_pos = light_pos*f_matrix;
  QVector4D ambient(0.4f, 0.4f, 0.4f, 0.4f);
  // Diffuse
  QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
  // Specular
  QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  vaos[Scene_edit_box_item_priv::Faces]->bind();
  d->program = &d->transparent_face_program;
  d->program->bind();
  d->program->setUniformValue("mvp_matrix", mvp_mat);
  d->program->setUniformValue("mv_matrix", mv_mat);
  d->program->setUniformValue("light_pos", light_pos);
  d->program->setUniformValue("light_diff",diffuse);
  d->program->setUniformValue("light_spec", specular);
  d->program->setUniformValue("light_amb", ambient);
  d->program->setUniformValue("spec_power", 51.8f);
  d->program->setAttributeValue("colors", QColor(128,128,128,128));
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->vertex_faces.size()/3));
  vaos[Scene_edit_box_item_priv::Faces]->release();
  d->program->release();
  glDisable(GL_BLEND);
  drawSpheres(viewer, f_matrix);
}

void Scene_edit_box_item::drawEdges(Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
  {
    d->computeElements();
    d->initializeBuffers(viewer);
  }
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }
  vaos[Scene_edit_box_item_priv::Edges]->bind();
  viewer->glLineWidth(6.0f);
  d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
  attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
  d->program->bind();
  d->program->setUniformValue("f_matrix", f_matrix);
  d->program->setAttributeValue("colors", QColor(Qt::black));
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->vertex_edges.size()/3));
  viewer->glLineWidth(1.0f);
  vaos[Scene_edit_box_item_priv::Edges]->release();
  d->program->release();
  if(renderingMode() == Wireframe)
  {
    drawSpheres(viewer, f_matrix);
  }
  drawHl(viewer);
}

void Scene_edit_box_item::compute_bbox() const
{
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();


  QVector3D min(d->pool[0], d->pool[1], d->pool[2]);
  QVector3D max(d->pool[3], d->pool[4], d->pool[5]);

  for(int i=0; i< 3; ++i)
  {
    min[i] += d->frame->translation()[i]-d->center_[i]-offset[i];
    max[i] += d->frame->translation()[i]-d->center_[i]-offset[i];
  }

  _bbox = Scene_item::Bbox(min.x(),min.y(),min.z(),max.x(),max.y(),max.z());
}


void
Scene_edit_box_item_priv::initializeBuffers(Viewer_interface *viewer)const
{

  //vao containing the data for the lines
  {
    program = item->getShaderProgram(Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT, viewer);
    program->bind();

    item->vaos[Edges]->bind();
    item->buffers[VertexEdges].bind();
    item->buffers[VertexEdges].allocate(vertex_edges.data(),
                                        static_cast<GLsizei>(vertex_edges.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    item->buffers[VertexEdges].release();
    item->vaos[Edges]->release();

    item->vaos[P_Edges]->bind();
    item->buffers[VertexEdges].bind();
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    item->buffers[VertexEdges].release();
    item->buffers[ColorsEdges].bind();
    item->buffers[ColorsEdges].allocate(color_edges.data(),
                                        static_cast<GLsizei>(color_edges.size()*sizeof(float)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",GL_FLOAT,0,3);
    item->buffers[ColorsEdges].release();
    item->vaos[P_Edges]->release();
    program->release();

  }

  //vao containing the data for the spheres
  {
    program = item->getShaderProgram(Scene_edit_box_item::PROGRAM_SPHERES, viewer);
    item->attribBuffers(viewer, Scene_edit_box_item::PROGRAM_SPHERES);

    program->bind();
    item->vaos[Spheres]->bind();
    item->buffers[VertexSpheres].bind();
    item->buffers[VertexSpheres].allocate(vertex_spheres.data(),
                                          static_cast<int>(vertex_spheres.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
    item->buffers[VertexSpheres].release();

    item->buffers[NormalSpheres].bind();
    item->buffers[NormalSpheres].allocate(normal_spheres.data(),
                                          static_cast<int>(normal_spheres.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
    item->buffers[NormalSpheres].release();

    item->buffers[CenterSpheres].bind();
    item->buffers[CenterSpheres].allocate(center_spheres.data(),
                                          static_cast<int>(center_spheres.size()*sizeof(float)));
    program->enableAttributeArray("center");
    program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
    item->buffers[CenterSpheres].release();
    program->disableAttributeArray("radius");
    program->disableAttributeArray("colors");

    viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
    viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);
    item->vaos[Spheres]->release();
    program->release();

    pick_sphere_program.bind();
    pick_sphere_program.bind();
    item->vaos[P_Spheres]->bind();
    item->buffers[VertexSpheres].bind();
    pick_sphere_program.enableAttributeArray("vertex");
    pick_sphere_program.setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
    item->buffers[VertexSpheres].release();

    item->buffers[NormalSpheres].bind();
    pick_sphere_program.disableAttributeArray("normals");
    pick_sphere_program.setAttributeValue("normals",QVector3D(0,0,0));
    item->buffers[NormalSpheres].release();

    item->buffers[CenterSpheres].bind();
    pick_sphere_program.enableAttributeArray("center");
    pick_sphere_program.setAttributeBuffer("center", GL_FLOAT, 0, 3);
    item->buffers[CenterSpheres].release();

    pick_sphere_program.disableAttributeArray("radius");
    viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);

    item->buffers[ColorsSpheres].bind();
    item->buffers[ColorsSpheres].allocate(color_spheres.data(),
                                          static_cast<int>(color_spheres.size()*sizeof(float)));
    pick_sphere_program.enableAttributeArray("colors");
    pick_sphere_program.setAttributeBuffer("colors", GL_FLOAT, 0, 3);
    item->buffers[ColorsSpheres].release();

    viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
    viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
    item->vaos[P_Spheres]->release();
    pick_sphere_program.release();
  }

  //vao containing the data for the faces
  {
    program = item->getShaderProgram(Scene_edit_box_item::PROGRAM_WITH_LIGHT, viewer);
    item->attribBuffers(viewer, Scene_edit_box_item::PROGRAM_WITH_LIGHT);

    program->bind();
    item->vaos[Faces]->bind();
    item->buffers[VertexFaces].bind();
    item->buffers[VertexFaces].allocate(vertex_faces.data(),
                                        static_cast<int>(vertex_faces.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
    item->buffers[VertexFaces].release();

    item->buffers[NormalFaces].bind();
    item->buffers[NormalFaces].allocate(normal_faces.data(),
                                        static_cast<int>(normal_faces.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
    item->buffers[NormalFaces].release();
    item->vaos[Faces]->release();
    program->release();


    program = item->getShaderProgram(Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT, viewer);
    item->attribBuffers(viewer, Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT);
    program->bind();
    item->vaos[P_Faces]->bind();
    item->buffers[VertexFaces].bind();
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
    item->buffers[VertexFaces].release();

    item->buffers[ColorsFaces].bind();
    item->buffers[ColorsFaces].allocate(color_faces.data(),
                                        static_cast<int>(color_faces.size()*sizeof(float)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors", GL_FLOAT, 0, 3);
    item->buffers[ColorsFaces].release();

    item->vaos[P_Faces]->release();
    program->release();
  }
  item->are_buffers_filled = true;
}

void push_xyz(std::vector<float> &v,
              const Scene_edit_box_item::Kernel::Point_3& p,
              qglviewer::Vec center_ = qglviewer::Vec(0,0,0))
{
  v.push_back(p.x()-center_.x);
  v.push_back(p.y()-center_.y);
  v.push_back(p.z()-center_.z);
}

void push_normal(std::vector<float> &v, int id)
{
  switch(id)
  {
  case 0:
    v.push_back(0);
    v.push_back(-1);
    v.push_back(0);
    break;
  case 1:
    v.push_back(-1);
    v.push_back(0);
    v.push_back(0);
    break;
  case 2:
    v.push_back(0);
    v.push_back(0);
    v.push_back(-1);
    break;
  case 3:
    v.push_back(1);
    v.push_back(0);
    v.push_back(0);
    break;
  case 4:
    v.push_back(0);
    v.push_back(0);
    v.push_back(1);
    break;
  case 5:
    v.push_back(0);
    v.push_back(1);
    v.push_back(0);
    break;
  default:
    break;
  }
}

void Scene_edit_box_item_priv::computeElements() const
{
  vertex_edges.clear();
  vertex_faces.clear();
  normal_faces.clear();
  center_spheres.clear();
  color_edges.clear();
  color_faces.clear();
  color_spheres.clear();

  //edges
  for( int i=0; i<12; ++i)
  {
    if(i<4)
    {
      vertex_edges.push_back(vertices[i].position().x()-center_.x);
      vertex_edges.push_back(vertices[i].position().y()-center_.y);
      vertex_edges.push_back(vertices[i].position().z()-center_.z);

      center_spheres.push_back(vertices[i].position().x()-center_.x);
      center_spheres.push_back(vertices[i].position().y()-center_.y);
      center_spheres.push_back(vertices[i].position().z()-center_.z);

      vertex_edges.push_back(vertices[i+4].position().x()-center_.x);
      vertex_edges.push_back(vertices[i+4].position().y()-center_.y);
      vertex_edges.push_back(vertices[i+4].position().z()-center_.z);

    }
    else if(i<8)
    {
      vertex_edges.push_back(vertices[i].position().x()-center_.x);
      vertex_edges.push_back(vertices[i].position().y()-center_.y);
      vertex_edges.push_back(vertices[i].position().z()-center_.z);

      center_spheres.push_back(vertices[i].position().x()-center_.x);
      center_spheres.push_back(vertices[i].position().y()-center_.y);
      center_spheres.push_back(vertices[i].position().z()-center_.z);

      vertex_edges.push_back(vertices[(i+1)%4 +4].position().x()-center_.x);
      vertex_edges.push_back(vertices[(i+1)%4 +4].position().y()-center_.y);
      vertex_edges.push_back(vertices[(i+1)%4 +4].position().z()-center_.z);
    }
    else
    {
      vertex_edges.push_back(vertices[i%4].position().x()-center_.x);
      vertex_edges.push_back(vertices[i%4].position().y()-center_.y);
      vertex_edges.push_back(vertices[i%4].position().z()-center_.z);

      vertex_edges.push_back(vertices[(i+1) %4].position().x()-center_.x);
      vertex_edges.push_back(vertices[(i+1) %4].position().y()-center_.y);
      vertex_edges.push_back(vertices[(i+1) %4].position().z()-center_.z);
    }
    color_edges.push_back(0);
    color_edges.push_back((20.0*i+10)/255);
    color_edges.push_back(0);

    color_edges.push_back(0);
    color_edges.push_back((20.0*i+10)/255);
    color_edges.push_back(0);

  }
  //faces
  for( int i=0; i<6; ++i)
  {
    push_xyz(vertex_faces, faces[i].vertices[0]->position(), center_);
    push_xyz(vertex_faces, faces[i].vertices[3]->position(), center_);
    push_xyz(vertex_faces, faces[i].vertices[2]->position(), center_);

    push_xyz(vertex_faces, faces[i].vertices[0]->position(), center_);
    push_xyz(vertex_faces, faces[i].vertices[2]->position(), center_);
    push_xyz(vertex_faces, faces[i].vertices[1]->position(), center_);

    for( int j=0; j<6; ++j)
    {
      push_normal(normal_faces, i);

      color_faces.push_back(0);
      color_faces.push_back(0);
      color_faces.push_back((20.0*i+10)/255);
    }
  }

  //spheres
  for( int i=0; i<8; ++i)
  {
    color_spheres.push_back((20.0*i+10)/255);
    color_spheres.push_back(0);
    color_spheres.push_back(0);
  }
}

Scene_edit_box_item::~Scene_edit_box_item()
{
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->setMouseTracking(false);

  delete d;
}

// Indicate if rendering mode is supported
bool Scene_edit_box_item::supportsRenderingMode(RenderingMode m) const {
  return (m==Wireframe || m==FlatPlusEdges);
}

Scene_item::ManipulatedFrame*
Scene_edit_box_item::manipulatedFrame()
{
  return d->frame;
}

double Scene_edit_box_item::point(short i, short j) const
{
  qglviewer::Vec pos(d->vertices[i].position().x()-d->center_.x,
                     d->vertices[i].position().y()-d->center_.y,
                     d->vertices[i].position().z()-d->center_.z);
  return (d->frame->inverseCoordinatesOf(pos))[j];
}

void Scene_edit_box_item::highlight()
{
  d->ready_to_hl = true;
  Viewer_interface* viewer = dynamic_cast<Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
  int type = -1, id = -1;
  //pick
  if(!d->selection_on)
  {
    d->picking(type, id, viewer);
    d->last_picked_id = id;
    d->last_picked_type = type;
  }
  //highlight
  d->hl_normal.clear();
  d->hl_vertex.clear();
  if(type !=-1)
  {
    switch(d->last_picked_type)
    {
    case 0:
    {
      //compute
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].position().x()-d->center_.x);
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].position().y()-d->center_.y);
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].position().z()-d->center_.z);
      //fill buffers
      d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_SPHERES, viewer);
      d->program->bind();

      vaos[Scene_edit_box_item_priv::S_Spheres]->bind();
      buffers[Scene_edit_box_item_priv::VertexSpheres].bind();
      buffers[Scene_edit_box_item_priv::VertexSpheres].allocate(d->vertex_spheres.data(),
                                      static_cast<int>(d->vertex_spheres.size()*sizeof(float)));
      d->program->enableAttributeArray("vertex");
      d->program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
      buffers[Scene_edit_box_item_priv::VertexSpheres].release();

      buffers[Scene_edit_box_item_priv::NormalSpheres].bind();
      buffers[Scene_edit_box_item_priv::NormalSpheres].allocate(d->normal_spheres.data(),
                                      static_cast<int>(d->normal_spheres.size()*sizeof(float)));
      d->program->enableAttributeArray("normals");
      d->program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
      buffers[Scene_edit_box_item_priv::NormalSpheres].release();

      buffers[Scene_edit_box_item_priv::S_Vertex].bind();
      buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                 static_cast<int>(d->hl_vertex.size()*sizeof(float)));
      d->program->enableAttributeArray("center");
      d->program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
      buffers[Scene_edit_box_item_priv::S_Vertex].release();
      d->program->disableAttributeArray("colors");

      viewer->glVertexAttribDivisor(d->program->attributeLocation("center"), 1);
      d->program->release();
      vaos[Scene_edit_box_item_priv::S_Spheres]->release();
      //draw
      d->hl_type = Scene_edit_box_item_priv::VERTEX;
      break;
    }
    case 1:
    {
      //compute
      d->hl_vertex.push_back(d->edges[d->last_picked_id].source->position().x()-d->center_.x);
      d->hl_vertex.push_back(d->edges[d->last_picked_id].source->position().y()-d->center_.y);
      d->hl_vertex.push_back(d->edges[d->last_picked_id].source->position().z()-d->center_.z);

      d->hl_vertex.push_back(d->edges[d->last_picked_id].target->position().x()-d->center_.x);
      d->hl_vertex.push_back(d->edges[d->last_picked_id].target->position().y()-d->center_.y);
      d->hl_vertex.push_back(d->edges[d->last_picked_id].target->position().z()-d->center_.z);

      //fill buffers
      d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT, viewer);
      d->program->bind();

      vaos[Scene_edit_box_item_priv::S_Edges]->bind();
      buffers[Scene_edit_box_item_priv::S_Vertex].bind();
      buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                 static_cast<GLsizei>(d->hl_vertex.size()*sizeof(float)));
      d->program->enableAttributeArray("vertex");
      d->program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
      buffers[Scene_edit_box_item_priv::S_Vertex].release();
      vaos[Scene_edit_box_item_priv::S_Edges]->release();
      d->program->release();
      //draw
      d->hl_type = Scene_edit_box_item_priv::EDGE;
      break;
    }
    case 2:
    {
      //compute
      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[0]->position(), d->center_);
      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[3]->position(), d->center_);
      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[2]->position(), d->center_);

      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[0]->position(), d->center_);
      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[2]->position(), d->center_);
      push_xyz(d->hl_vertex, d->faces[d->last_picked_id].vertices[1]->position(), d->center_);

      for( int j=0; j<6; ++j)
      {
        push_normal(d->hl_normal, d->last_picked_id);
      }
      //fill buffers
      d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_WITH_LIGHT, viewer);
      attribBuffers(viewer, Scene_edit_box_item::PROGRAM_WITH_LIGHT);

      d->program->bind();
      vaos[Scene_edit_box_item_priv::S_Faces]->bind();
      buffers[Scene_edit_box_item_priv::S_Vertex].bind();
      buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                 static_cast<int>(d->hl_vertex.size()*sizeof(float)));
      d->program->enableAttributeArray("vertex");
      d->program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
      buffers[Scene_edit_box_item_priv::S_Normal].release();

      buffers[Scene_edit_box_item_priv::S_Normal].bind();
      buffers[Scene_edit_box_item_priv::S_Normal].allocate(d->hl_normal.data(),
                                 static_cast<int>(d->hl_normal.size()*sizeof(float)));
      d->program->enableAttributeArray("normals");
      d->program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
      buffers[Scene_edit_box_item_priv::S_Normal].release();
      vaos[Scene_edit_box_item_priv::S_Faces]->release();
      d->program->release();

      //draw
      d->hl_type = Scene_edit_box_item_priv::FACE;
      break;
    }
    default:
      d->hl_type = Scene_edit_box_item_priv::NO_TYPE;
      break;
    }
  }
  itemChanged();

  d->ready_to_hl = false;
}

void Scene_edit_box_item::clearHL()
{
  Viewer_interface* viewer = dynamic_cast<Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
  d->hl_normal.clear();
  d->hl_vertex.clear();

  d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_SPHERES, viewer);
  d->program->bind();

  vaos[Scene_edit_box_item_priv::S_Spheres]->bind();
  buffers[Scene_edit_box_item_priv::VertexSpheres].bind();
  buffers[Scene_edit_box_item_priv::VertexSpheres].allocate(d->vertex_spheres.data(),
                                                            static_cast<int>(d->vertex_spheres.size()*sizeof(float)));
  d->program->enableAttributeArray("vertex");
  d->program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  buffers[Scene_edit_box_item_priv::VertexSpheres].release();

  buffers[Scene_edit_box_item_priv::NormalSpheres].bind();
  buffers[Scene_edit_box_item_priv::NormalSpheres].allocate(d->normal_spheres.data(),
                                                            static_cast<int>(d->normal_spheres.size()*sizeof(float)));
  d->program->enableAttributeArray("normals");
  d->program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  buffers[Scene_edit_box_item_priv::NormalSpheres].release();

  buffers[Scene_edit_box_item_priv::S_Vertex].bind();
  buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                                       static_cast<int>(d->hl_vertex.size()*sizeof(float)));
  d->program->enableAttributeArray("center");
  d->program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
  buffers[Scene_edit_box_item_priv::S_Vertex].release();
  d->program->disableAttributeArray("colors");

  viewer->glVertexAttribDivisor(d->program->attributeLocation("center"), 1);
  d->program->release();
  vaos[Scene_edit_box_item_priv::S_Spheres]->release();
  //draw
  d->hl_type = Scene_edit_box_item_priv::VERTEX;
  d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT, viewer);
  d->program->bind();

  vaos[Scene_edit_box_item_priv::S_Edges]->bind();
  buffers[Scene_edit_box_item_priv::S_Vertex].bind();
  buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                                       static_cast<GLsizei>(d->hl_vertex.size()*sizeof(float)));
  d->program->enableAttributeArray("vertex");
  d->program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  buffers[Scene_edit_box_item_priv::S_Vertex].release();
  vaos[Scene_edit_box_item_priv::S_Edges]->release();
  d->program->release();

  d->program = getShaderProgram(Scene_edit_box_item::PROGRAM_WITH_LIGHT, viewer);
  attribBuffers(viewer, Scene_edit_box_item::PROGRAM_WITH_LIGHT);

  d->program->bind();
  vaos[Scene_edit_box_item_priv::S_Faces]->bind();
  buffers[Scene_edit_box_item_priv::S_Vertex].bind();
  buffers[Scene_edit_box_item_priv::S_Vertex].allocate(d->hl_vertex.data(),
                                                       static_cast<int>(d->hl_vertex.size()*sizeof(float)));
  d->program->enableAttributeArray("vertex");
  d->program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  buffers[Scene_edit_box_item_priv::S_Normal].release();

  buffers[Scene_edit_box_item_priv::S_Normal].bind();
  buffers[Scene_edit_box_item_priv::S_Normal].allocate(d->hl_normal.data(),
                                                       static_cast<int>(d->hl_normal.size()*sizeof(float)));
  d->program->enableAttributeArray("normals");
  d->program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  buffers[Scene_edit_box_item_priv::S_Normal].release();
  vaos[Scene_edit_box_item_priv::S_Faces]->release();
  d->program->release();
  d->hl_type = Scene_edit_box_item_priv::NO_TYPE;

  itemChanged();

}
void Scene_edit_box_item_priv::reset_selection()
{
  selection_on = false;
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->setManipulatedFrame(frame);
  viewer->setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, QGLViewer::SELECT);
  constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
  selected_vertices.clear();
}

//intercept events for picking
bool Scene_edit_box_item::eventFilter(QObject *, QEvent *event)
{
  if(event->type() == QEvent::MouseButtonPress)
  {
    QMouseEvent* e = static_cast<QMouseEvent*>(event);
    if(e->modifiers() == Qt::ShiftModifier)
    {
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      Viewer_interface* v_i = dynamic_cast<Viewer_interface*>(viewer);
      //pick
      int type, picked;
      d->picked_pixel = e->pos();
      d->picking(type, picked, v_i);
      if(type !=-1)
      {
        bool found = false;
        QApplication::setOverrideCursor(Qt::DragMoveCursor);
        qglviewer::Vec pos = viewer->camera()->pointUnderPixel(d->picked_pixel, found);
        if(found)
        {
          d->rf_last_pos = pos;
          d->remodel_frame->setPosition(pos);
        }
        for(int i=0; i<8; ++i)
          for(int j=0; j<3; ++j)
            d->last_pool[i][j] = d->vertices[i][j];
        d->selection_on = true;
        if(type == 0)
        {
          d->selected_vertices.push_back(&d->vertices[picked]);
          d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
          d->remodel_frame->setConstraint(&d->constraint);
        }
        else if(type == 1)
        {
          d->selected_vertices.push_back(d->edges[picked].source);
          d->selected_vertices.push_back(d->edges[picked].target);
          Kernel::Point_3 s(d->edges[picked].source->position()), t(d->edges[picked].target->position());

          qglviewer::Vec normal(t.x()-s.x(), t.y()-s.y(), t.z()-s.z());
          d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::PLANE);
          d->constraint.setTranslationConstraintDirection(normal);
          d->remodel_frame->setConstraint(&d->constraint);
        }
        else if(type == 2)
        {
          for(int i=0; i<4; ++i)
            d->selected_vertices.push_back(d->faces[picked].vertices[i]);
          Kernel::Point_3 a1(d->faces[picked].vertices[1]->position()), a0(d->faces[picked].vertices[0]->position())
              ,a3(d->faces[picked].vertices[3]->position());
          QVector3D a(a1.x()-a0.x(), a1.y()-a0.y(),a1.z()-a0.z()),b(a3.x()-a0.x(), a3.y()-a0.y(),a3.z()-a0.z());
          QVector3D n = QVector3D::crossProduct(a,b);

          d->remodel_frame->setConstraint(&d->constraint);
          d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
          d->constraint.setTranslationConstraintDirection(qglviewer::Vec(n.x(), n.y(), n.z()));
        }

        viewer->setManipulatedFrame(d->remodel_frame);
        viewer->setMouseBinding(
              Qt::ShiftModifier,
              Qt::LeftButton,
              QGLViewer::FRAME,
              QGLViewer::TRANSLATE);
      }
      else
      {
        d->reset_selection();
      }
    }
    return false;
  }
  else if(event->type() == QEvent::MouseMove)
  {
    QMouseEvent* e = static_cast<QMouseEvent*>(event);
    if(e->modifiers() == Qt::ShiftModifier)
    {
      if(d->selection_on)
      {
        d->remodel_frame->setOrientation(d->frame->orientation());
        qglviewer::Vec td(d->remodel_frame->transformOf(d->remodel_frame->position() - d->rf_last_pos));
        QVector3D dir(td.x, td.y, td.z);
        d->remodel_box(dir);
      }
      d->ready_to_hl= true;
      d->picked_pixel = e->pos();
      QTimer::singleShot(0, this, SLOT(highlight()));
    }
    else if(e->modifiers() == Qt::ControlModifier &&
            e->buttons() == Qt::LeftButton)
    {
      QApplication::setOverrideCursor(d->rotate_cursor);
    }
    else if(d->selection_on)
    {
      d->reset_selection();
    }
    d->picked_pixel = e->pos();
    return false;
  }
  else if(event->type() == QEvent::MouseButtonRelease)
  {
    d->reset_selection();
    QApplication::setOverrideCursor(QCursor());
  }

  else if(event->type() == QEvent::KeyPress)
  {
     QKeyEvent* e = static_cast<QKeyEvent*>(event);
     if(e->key() == Qt::Key_Shift)
     {
       d->ready_to_hl= true;
       QTimer::singleShot(0, this, SLOT(highlight()));
     }
  }
  else if(event->type() == QEvent::KeyRelease)
  {
     QKeyEvent* e = static_cast<QKeyEvent*>(event);
     if(e->key() == Qt::Key_Shift)
       QTimer::singleShot(0, this, SLOT(clearHL()));
     else if(e->key() == Qt::Key_Control)
     {

       QApplication::setOverrideCursor(QCursor());
     }
  }
  return false;
}

void Scene_edit_box_item_priv::draw_picking(Viewer_interface* viewer)
{

  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)frame->matrix()[i];
  }
  GLdouble d_mat[16];
  QMatrix4x4 mvp_mat;
  viewer->camera()->getModelViewProjectionMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mvp_mat.data()[i] = GLfloat(d_mat[i]);
  mvp_mat = mvp_mat*f_matrix;
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  QVector4D light_pos(0.0f,0.0f,1.0f, 1.0f );
  light_pos = light_pos*f_matrix;

  if(item->renderingMode() == FlatPlusEdges)
  {
    item->vaos[P_Faces]->bind();
    program = item->getShaderProgram(Scene_item::PROGRAM_WITHOUT_LIGHT, viewer);
    item->attribBuffers(viewer, Scene_item::PROGRAM_WITHOUT_LIGHT);
    program->bind();
    program->setUniformValue("mvp_matrix", mvp_mat);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertex_faces.size()/3));
    item->vaos[P_Faces]->release();
    program->release();
  }
  item->vaos[P_Spheres]->bind();
  pick_sphere_program.bind();
  pick_sphere_program.setUniformValue("mvp_matrix", mvp_mat);
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(vertex_spheres.size()/3),
                                static_cast<GLsizei>(8));
  pick_sphere_program.release();
  item->vaos[P_Spheres]->release();

  item->vaos[P_Edges]->bind();
  viewer->glLineWidth(6.0f);
  program = item->getShaderProgram(Scene_item::PROGRAM_WITHOUT_LIGHT);
  item->attribBuffers(viewer, Scene_item::PROGRAM_WITHOUT_LIGHT);
  program->bind();
  program->setUniformValue("f_matrix", f_matrix);
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(vertex_edges.size()/3));
  viewer->glLineWidth(1.0f);
  item->vaos[P_Edges]->release();
  program->release();
}

void Scene_edit_box_item_priv::remodel_box(const QVector3D &dir)
{
  qglviewer::AxisPlaneConstraint::Type prev_cons = constraint.translationConstraintType();
  constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
  Q_FOREACH(Scene_edit_box_item::vertex*  selected_vertex, selected_vertices )
  {
    int id = selected_vertex->id;
    *selected_vertex->x = applyX(id, last_pool[id][0], dir.x());
    *selected_vertex->y = applyY(id, last_pool[id][1], dir.y());
    *selected_vertex->z = applyZ(id, last_pool[id][2], dir.z());
    for( int i=0; i<3; ++i)
      relative_center_[i] =(pool[i]+pool[i+3])/2 - center_[i];
    for( int i=0; i<3; ++i)
      center_[i] =(pool[i]+pool[i+3])/2;
    frame->translate(frame->inverseTransformOf(relative_center_));
  }
  item->invalidateOpenGLBuffers();
  constraint.setTranslationConstraintType(prev_cons);
}

double Scene_edit_box_item_priv::applyX(int id, double x, double dirx)
{
  switch(id)
  {
  case 0:
  case 1:
  case 4:
  case 5:
    if(x+dirx < pool[3])
      return x+dirx;
    else
      return pool[3];
  case 2:
  case 3:
  case 6:
  case 7:
    if(x+dirx > pool[0])
      return x+dirx;
    else
      return pool[0];
  default:
    return 0;
  }
  return 0;
}

double Scene_edit_box_item_priv::applyY(int id, double y, double diry)
{
  switch(id)
  {
  case 0:
  case 1:
  case 2:
  case 3:
    if(y+diry < pool[4])
      return y+diry;
    else
      return pool[4];
  case 4:
  case 5:
  case 6:
  case 7:
    if(y+diry > pool[1])
      return y+diry;
    else
      return pool[1];
  default:
    return 0;
  }
  return 0;
}

double Scene_edit_box_item_priv::applyZ(int id, double z, double dirz)
{
  switch(id)
  {
  case 1:
  case 2:
  case 5:
  case 6:
    if(z+dirz < pool[5])
      return z+dirz;
    else
      return pool[5];
  case 0:
  case 3:
  case 4:
  case 7:
    if(z+dirz > pool[2])
      return z+dirz;
    else
      return pool[2];
  default:
    return 0;
  }
  return 0;
}

//type : 0 = vertex, 1 = edge, 2 = face
void Scene_edit_box_item_priv::picking(int& type, int& id, Viewer_interface *viewer)
{
  type = -1;
  id = -1;
  int deviceWidth = viewer->camera()->screenWidth();
  int deviceHeight = viewer->camera()->screenHeight();
  QOpenGLFramebufferObject* fbo = new QOpenGLFramebufferObject(deviceWidth, deviceHeight,QOpenGLFramebufferObject::Depth);
  fbo->bind();
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  QColor bgColor(viewer->backgroundColor());
  //draws the image in the fbo
  viewer->setBackgroundColor(::Qt::white);
  draw_picking(viewer);

  int rowLength = deviceWidth * 4; // data asked in RGBA,so 4 bytes.
  const static int dataLength = rowLength * deviceHeight;
  GLubyte* buffer = new GLubyte[dataLength];
  // Qt uses upper corner for its origin while GL uses the lower corner.
  glReadPixels(picked_pixel.x(), deviceHeight-1-picked_pixel.y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
  //decode ID and pick (don't forget the case nothing is picked
  if(!(buffer[0]==buffer[1] && buffer[1]==buffer[2]))
  {
    int r(std::ceil((buffer[0]-10)/20)), g(std::ceil((buffer[1]-10)/20)), b(std::ceil((buffer[2]-10)/20));
    id = (std::max)(r,g);
    id = (std::max)(id,b);
    if(buffer[0] > 0)
    {
      if(id <8)
        type = 0 ;
    }
    else if(buffer[1] > 0)
    {
      if(id <12)
      {
        type = 1;
      }
    }
    else if(buffer[2] > 0)
    {
      if(id <6)
      {
        type = 2;
      }
    }
  }
  delete[] buffer;
  viewer->setBackgroundColor(bgColor);
  fbo->release();
  delete fbo;
}

void Scene_edit_box_item::drawHl(Viewer_interface* viewer)const
{
  GLfloat offset_factor;
  GLfloat offset_units;
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }
  GLdouble d_mat[16];
  QMatrix4x4 mvp_mat;
  viewer->camera()->getModelViewProjectionMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mvp_mat.data()[i] = GLfloat(d_mat[i]);
  mvp_mat = mvp_mat*f_matrix;
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  QVector4D light_pos(0.0f,0.0f,1.0f, 1.0f );
  light_pos = light_pos*f_matrix;
  QVector4D ambient(0.4f, 0.4f, 0.4f, 0.4f);
  // Diffuse
  QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
  // Specular
  QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);

  if(d->hl_type == Scene_edit_box_item_priv::VERTEX)
  {
    vaos[Scene_edit_box_item_priv::S_Spheres]->bind();
    attribBuffers(viewer, PROGRAM_SPHERES);
    d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
    d->program->bind();
    d->program->setUniformValue("mvp_matrix", mvp_mat);
    d->program->setUniformValue("mv_matrix", mv_mat);
    d->program->setUniformValue("light_pos", light_pos);
    d->program->setAttributeValue("colors", QColor(Qt::yellow));
    double radius =std::sqrt(
        (point(6,0) - point(0,0)) * (point(6,0) - point(0,0)) +
        (point(6,1) - point(0,1)) * (point(6,1) - point(0,1)) +
        (point(6,2) - point(0,2)) * (point(6,2) - point(0,2))) *0.02 ;
    d->program->setAttributeValue("radius", radius);
    viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                  static_cast<GLsizei>(d->vertex_spheres.size()/3),
                                  static_cast<GLsizei>(d->hl_vertex.size()/3));

    d->program->release();
    vaos[Scene_edit_box_item_priv::S_Spheres]->release();
  }
  else if(d->hl_type == Scene_edit_box_item_priv::EDGE)
  {
    vaos[Scene_edit_box_item_priv::S_Edges]->bind();
    viewer->glLineWidth(6.0f);
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    d->program->setUniformValue("f_matrix", f_matrix);
    d->program->setAttributeValue("colors", QColor(Qt::yellow));
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->hl_vertex.size()/3));
    viewer->glLineWidth(1.0f);
    vaos[Scene_edit_box_item_priv::S_Edges]->release();
    d->program->release();
  }
  else if(d->hl_type == Scene_edit_box_item_priv::FACE)
  {
    viewer->glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    vaos[Scene_edit_box_item_priv::S_Faces]->bind();
    d->program = &d->transparent_face_program;
    d->program->bind();
    d->program->setUniformValue("mvp_matrix", mvp_mat);
    d->program->setUniformValue("mv_matrix", mv_mat);
    d->program->setUniformValue("light_pos", light_pos);
    d->program->setUniformValue("light_diff",diffuse);
    d->program->setUniformValue("light_spec", specular);
    d->program->setUniformValue("light_amb", ambient);
    d->program->setUniformValue("spec_power", 51.8f);
    d->program->setAttributeValue("colors", QColor(128,128,0,128));
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->hl_vertex.size()/3));
    vaos[Scene_edit_box_item_priv::S_Faces]->release();
    d->program->release();
    glPolygonOffset(offset_factor, offset_units);
    glDisable(GL_BLEND);

  }
}
