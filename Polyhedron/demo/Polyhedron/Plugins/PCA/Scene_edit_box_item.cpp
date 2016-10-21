#include "Scene_edit_box_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLShaderProgram>
using namespace CGAL::Three;
struct Scene_edit_box_item::vertex{
  short id;
  face* face_;
  double *x;
  double *y;
  double *z;

  CGAL::Point_3<Kernel> position()const
  {
   return  CGAL::Point_3<Kernel>(*x,*y,*z);
  }
};
struct Scene_edit_box_item::edge{
  short id;
  vertex* source;
  vertex* target;
  face* face_;
};
struct Scene_edit_box_item::face{
  short id;
  edge* edges[4];
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
    Arrow,
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
    VertexArrow,
    NormalArrow,
    NumberOfVbos
  };

  Scene_edit_box_item_priv(const Scene_interface *scene_interface, Scene_edit_box_item* ebi)
  {
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
    frame->setPosition(center_);
    frame->setSpinningSensitivity(100.0); //forbid spinning
    constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    constraint.setRotationConstraintDirection(qglviewer::Vec(.0,.0,.1));
    frame->setConstraint(&constraint);
    //create the sphere model
    pool[0] = bb.xmin(); pool[3] = bb.xmax();
    pool[1] = bb.ymin(); pool[4] = bb.ymax();
    pool[2] = bb.zmin(); pool[5] = bb.zmax();

    double sq_diag =
        (bb.xmax() - bb.xmin()) * (bb.xmax() - bb.xmin()) +
        (bb.ymax() - bb.ymin()) * (bb.ymax() - bb.ymin()) +
        (bb.zmax() - bb.zmin()) * (bb.zmax() - bb.zmin()) ;

    vertex_spheres.resize(0);
    normal_spheres.resize(0);
    create_flat_sphere((float)(0.01*sqrt(sq_diag)), vertex_spheres, normal_spheres,10);
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
      edges[i].id = i;
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
      faces[0].edges[(i+1)%4] = &edges[i+7];
    }
    faces[0].id =0;

    for( int i=1; i<4; ++i)
    {
      faces[i].vertices[0] = &vertices[i];
      faces[i].vertices[1] = &vertices[(i-1)];
      faces[i].vertices[2] = &vertices[i+3];
      faces[i].vertices[3] = &vertices[i+4];

      faces[i].edges[0] = &edges[i+7];
      faces[i].edges[1] = &edges[i-1];
      faces[i].edges[2] = &edges[i+3];
      faces[i].edges[3] = &edges[i];
      faces[i].id = i;
    }

    faces[4].vertices[0] = &vertices[0];
    faces[4].vertices[1] = &vertices[3];
    faces[4].vertices[2] = &vertices[7];
    faces[4].vertices[3] = &vertices[4];

    faces[4].edges[0] = &edges[0];
    faces[4].edges[1] = &edges[3];
    faces[4].edges[2] = &edges[7];
    faces[4].edges[3] = &edges[4];
    faces[4].id =4;

    for( int i=0; i<4; ++i)
    {
      faces[5].vertices[i] = &vertices[i+4];
      faces[5].edges[i] = &edges[i+4];
    }
    faces[5].id =5;

    vertex_faces.resize(0);
    normal_faces.resize(0);

    for( int i=0; i<6; ++i)
    {
      for( int j=0; j<4; ++j)
      {
        faces[i].vertices[j]->face_ = &faces[i];
        faces[i].edges[j]->face_ = &faces[i];
      }
    }

    pick_sphere_program.addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_spheres.v");
    pick_sphere_program.addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_without_light.f");
    pick_sphere_program.bindAttributeLocation("colors", 1);
    pick_sphere_program.link();
  }

  mutable std::vector<float> vertex_edges;
  mutable std::vector<float> color_edges;
  mutable std::vector<float> vertex_spheres;
  mutable std::vector<float> normal_spheres;
  mutable std::vector<float> center_spheres;
  mutable std::vector<float> color_spheres;
  mutable std::vector<float> vertex_arrow;
  mutable std::vector<float> normal_arrow;
  mutable std::vector<float> vertex_faces;
  mutable std::vector<float> normal_faces;
  mutable std::vector<float> color_faces;
  double pool[6];


  qglviewer::ManipulatedFrame* frame;
  qglviewer::ManipulatedFrame* remodel_frame;
  qglviewer::Vec rf_last_pos;
  qglviewer::LocalConstraint constraint;
  qglviewer::Vec center_;
  qglviewer::Vec relative_center_;

  mutable QOpenGLShaderProgram pick_sphere_program;
  mutable Scene_edit_box_item::vertex vertices[8];
  mutable Scene_edit_box_item::edge edges[12];
  mutable Scene_edit_box_item::face faces[6];

  Scene_edit_box_item::vertex* selected_vertex;
  Scene_edit_box_item::edge* selected_edge;
  Scene_edit_box_item::face* selected_face;
  void reset_selection();
  bool selection_on;

  mutable QOpenGLShaderProgram* program;
  void initializeBuffers(Viewer_interface *viewer)const;

  void computeElements() const;
  void draw_picking(Viewer_interface*);
  void remodel_box(const QVector3D &dir);

  bool applyX(double x, double dirx);
  bool applyY(double y, double diry);
  bool applyZ(double z, double dirz);
  const Scene_interface* scene;
  Scene_edit_box_item* item;
};


Scene_edit_box_item::Scene_edit_box_item(const Scene_interface *scene_interface)
  :  Scene_item(NumberOfVbos,NumberOfVaos)

{
  d = new Scene_edit_box_item_priv(scene_interface, this);

  are_buffers_filled = false;
}
QString Scene_edit_box_item::toolTip() const {

  return QString();
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


  vaos[Scene_edit_box_item_priv::Faces]->bind();
  d->program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
  attribBuffers(viewer, PROGRAM_WITH_LIGHT);
  d->program->bind();
  d->program->setUniformValue("mvp_matrix", mvp_mat);
  d->program->setUniformValue("mv_matrix", mv_mat);
  d->program->setUniformValue("light_pos", light_pos);
  d->program->setAttributeValue("colors", QColor(this->color()));
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->vertex_faces.size()/3));
  vaos[Scene_edit_box_item_priv::Faces]->release();
  d->program->release();

  vaos[Scene_edit_box_item_priv::Spheres]->bind();
  d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
  attribBuffers(viewer, PROGRAM_SPHERES);
  d->program->bind();
  d->program->setUniformValue("mvp_matrix", mvp_mat);
  d->program->setUniformValue("light_pos", light_pos);
  d->program->setAttributeValue("colors", QColor(Qt::red));
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(d->vertex_spheres.size()/3),
                                static_cast<GLsizei>(8));
  d->program->release();
  vaos[Scene_edit_box_item_priv::Spheres]->release();
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
  vaos[Edges]->bind();
  viewer->glLineWidth(4.0f);
  d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
  attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
  d->program->bind();
  d->program->setUniformValue("f_matrix", f_matrix);
  d->program->setAttributeValue("colors", QColor(Qt::black));
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->vertex_edges.size()/3));
  viewer->glLineWidth(1.0f);
  vaos[Edges]->release();
  d->program->release();
  if(renderingMode() == Wireframe)
  {
    if (!are_buffers_filled)
    {
      d->computeElements();
      d->initializeBuffers(viewer);
    }
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


    d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attribBuffers(viewer, PROGRAM_SPHERES);
    d->program->bind();
    d->program->setUniformValue("mvp_matrix", mvp_mat);
    d->program->setUniformValue("mv_matrix", mv_mat);
    d->program->setUniformValue("light_pos", light_pos);
    d->program->setAttributeValue("colors", QColor(Qt::red));
    viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                  static_cast<GLsizei>(d->vertex_spheres.size()/3),
                                  static_cast<GLsizei>(8));
    d->program->release();
    vaos[Scene_edit_box_item_priv::Spheres]->release();
  }
}

void Scene_edit_box_item::compute_bbox() const
{
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }

  double xmin,ymin,zmin,xmax,ymax,zmax;
  xmin =d->vertices[0].position().x();
  ymin =d->vertices[0].position().y();
  zmin =d->vertices[0].position().z();
  xmax =d->vertices[6].position().x();
  ymax =d->vertices[6].position().y();
  zmax =d->vertices[6].position().z();
  QVector3D min(xmin, ymin, zmin);
  QVector3D max(xmax, ymax, zmax);
  min = f_matrix*min;
  max = f_matrix*max;
  Scene_item::Bbox bb(min.x(),min.y(),min.z(),max.x(),max.y(),max.z());
  _bbox = bb;
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
    program->enableAttributeArray("center");
    item->buffers[CenterSpheres].allocate(center_spheres.data(),
                                          static_cast<int>(center_spheres.size()*sizeof(float)));
    program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
    item->buffers[VertexEdges].release();
    program->disableAttributeArray("radius");
    program->setAttributeValue("radius",1);
    program->disableAttributeArray("colors");

    viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
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
    pick_sphere_program.setAttributeValue("radius",1);

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
  QApplication::setOverrideCursor(Qt::WaitCursor);
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

  QApplication::restoreOverrideCursor();
}

Scene_edit_box_item::~Scene_edit_box_item()
{
  delete d;
}

// Indicate if rendering mode is supported
bool Scene_edit_box_item::supportsRenderingMode(RenderingMode m) const {
  switch(m)
  {
  case Wireframe:
  case FlatPlusEdges:
    return true;
  default:
    return false;
  }
}

Scene_item::ManipulatedFrame*
Scene_edit_box_item::manipulatedFrame()
{
  return d->frame;
}

double Scene_edit_box_item::point(short i, short j) const
{
  QVector3D pos(d->vertices[i].position().x()-d->center_.x,
                d->vertices[i].position().y()-d->center_.y,
                d->vertices[i].position().z()-d->center_.z);
  QMatrix4x4 f_matrix;
  for (int k=0; k<16; ++k){
    f_matrix.data()[k] = (float)d->frame->matrix()[k];
  }
  return (f_matrix*pos)[j];
}

void Scene_edit_box_item::highlight()
{
  /*  if(is_ready_to_highlight)
  {
    // highlight with mouse move event
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    bool found = false;
    const qglviewer::Vec& point = camera->pointUnderPixel(hl_pos, found);
    if(found)
    {
      const qglviewer::Vec& orig = camera->position();
      const qglviewer::Vec& dir = point - orig;
      is_highlighting = true;

      is_highlighting = false;
    }
    else
    {
      Q_EMIT clearHL();
    }
    is_ready_to_highlight = false;
  }*/
}

void Scene_edit_box_item_priv::reset_selection()
{
  selection_on = false;
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->setManipulatedFrame(frame);
  viewer->setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, QGLViewer::SELECT);
  constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
  selected_vertex= NULL;
  selected_edge= NULL;
  selected_face= NULL;
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
      int deviceWidth = viewer->camera()->screenWidth();
      int deviceHeight = viewer->camera()->screenHeight();
      QOpenGLFramebufferObject* fbo = new QOpenGLFramebufferObject(deviceWidth, deviceHeight,QOpenGLFramebufferObject::Depth);
      fbo->bind();
      glEnable(GL_DEPTH_TEST);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      QColor bgColor(viewer->backgroundColor());
      //draws the image in the fbo
      viewer->setBackgroundColor(::Qt::white);
      Viewer_interface* v_i = dynamic_cast<Viewer_interface*>(viewer);
      //draw_picking
      d->draw_picking(v_i);

      int rowLength = deviceWidth * 4; // data asked in RGBA,so 4 bytes.
      const static int dataLength = rowLength * deviceHeight;
      GLubyte* buffer = new GLubyte[dataLength];
      // Qt uses upper corner for its origin while GL uses the lower corner.
      glReadPixels(e->pos().x(), deviceHeight-1-e->pos().y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
      //decode ID and pick (don't forget the case nothing is picked
      if(!(buffer[0]==buffer[1] && buffer[1]==buffer[2]))
      {
        d->selection_on = true;
        d->rf_last_pos = d->remodel_frame->position();
        int r(std::ceil((buffer[0]-10)/20)), g(std::ceil((buffer[1]-10)/20)), b(std::ceil((buffer[2]-10)/20));
        int picked = (std::max)(r,g);
        picked = (std::max)(picked,b);
        if(buffer[0] > 0)
        {
          if(picked <8)
          {
            d->selected_vertex = &d->vertices[picked];
            d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
            d->remodel_frame->setConstraint(&d->constraint);
          }
        }
        else if(buffer[1] > 0)
        {
          if(picked <12)
          {
            d->selected_edge = &d->edges[picked];
            Kernel::Point_3 s(d->selected_edge->source->position()), t(d->selected_edge->target->position());

            qglviewer::Vec normal(t.x()-s.x(), t.y()-s.y(), t.z()-s.z());
            d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::PLANE);
            d->constraint.setTranslationConstraintDirection(normal);
            d->remodel_frame->setConstraint(&d->constraint);
          }
        }
        else if(buffer[2] > 0)
        {
          if(picked <6)
          {
            d->selected_face= &d->faces[picked];
            Kernel::Point_3 a1(d->selected_face->vertices[1]->position()), a0(d->selected_face->vertices[0]->position())
                ,a3(d->selected_face->vertices[3]->position());
            QVector3D a(a1.x()-a0.x(), a1.y()-a0.y(),a1.z()-a0.z()),b(a3.x()-a0.x(), a3.y()-a0.y(),a3.z()-a0.z());
            QVector3D n = QVector3D::crossProduct(a,b);

            d->remodel_frame->setConstraint(&d->constraint);
            d->constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
            d->constraint.setTranslationConstraintDirection(qglviewer::Vec(n.x(), n.y(), n.z()));
          }
        }
        viewer->setBackgroundColor(bgColor);
        fbo->release();
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
      delete fbo;
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
        d->rf_last_pos = d->remodel_frame->position();
        d->remodel_box(dir);
        return false;
      }
    }
    else if(d->selection_on)
    {
      d->reset_selection();
    }
  }
  else if(event->type() == QEvent::MouseButtonRelease)
  {
    d->reset_selection();
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
  viewer->glLineWidth(4.0f);
  program = item->getShaderProgram(Scene_item::PROGRAM_WITHOUT_LIGHT);
  item->attribBuffers(viewer, Scene_item::PROGRAM_WITHOUT_LIGHT);
  program->bind();
  program->setUniformValue("f_matrix", f_matrix);
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(vertex_edges.size()/3));
  viewer->glLineWidth(1.0f);
  item->vaos[P_Edges]->release();
  program->release();
}

//!\todo redo the API to only use a list of selected vertices
void Scene_edit_box_item_priv::remodel_box(const QVector3D &dir)
{
  qglviewer::AxisPlaneConstraint::Type prev_cons = constraint.translationConstraintType();
  constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);

  if(selected_vertex != NULL)
  {
    if(applyX(*selected_vertex->x, dir.x()))
       *selected_vertex->x += dir.x();
    if(applyY(*selected_vertex->y, dir.y()))
      *selected_vertex->y += dir.y();
    if(applyZ(*selected_vertex->z, dir.z()))
      *selected_vertex->z += dir.z();
    for( int i=0; i<3; ++i)
      relative_center_[i] =(pool[i]+pool[i+3])/2 - center_[i];
    for( int i=0; i<3; ++i)
      center_[i] =(pool[i]+pool[i+3])/2;
    frame->translate(frame->inverseTransformOf(relative_center_));
    item->invalidateOpenGLBuffers();
  }
  else if(selected_edge != NULL)
  {
    if(applyX(*selected_edge->source->x, dir.x()))
       *selected_edge->source->x += dir.x();
    if(applyY(*selected_edge->source->y, dir.y()))
      *selected_edge->source->y += dir.y();
    if(applyZ(*selected_edge->source->z, dir.z()))
      *selected_edge->source->z += dir.z();
    if(applyX(*selected_edge->target->x, dir.x()))
       *selected_edge->target->x += dir.x();
    if(applyY(*selected_edge->target->y, dir.y()))
      *selected_edge->target->y += dir.y();
    if(applyZ(*selected_edge->target->z, dir.z()))
      *selected_edge->target->z += dir.z();
    for( int i=0; i<3; ++i)
      relative_center_[i] =(pool[i]+pool[i+3])/2 - center_[i];
    for( int i=0; i<3; ++i)
      center_[i] =(pool[i]+pool[i+3])/2;
    frame->translate(frame->inverseTransformOf(relative_center_));
    item->invalidateOpenGLBuffers();
  }
  else if(selected_face != NULL)
  {
    for( int i=0; i<4; ++i)
    {
      if(applyX(*selected_face->vertices[i]->x, dir.x()))
         *selected_face->vertices[i]->x += dir.x();
      if(applyY(*selected_face->vertices[i]->y, dir.y()))
        *selected_face->vertices[i]->y += dir.y();
      if(applyZ(*selected_face->vertices[i]->z, dir.z()))
        *selected_face->vertices[i]->z += dir.z();
    }
    for( int i=0; i<3; ++i)
      relative_center_[i] =(pool[i]+pool[i+3])/2 - center_[i];
    for( int i=0; i<3; ++i)
      center_[i] =(pool[i]+pool[i+3])/2;
    frame->translate(frame->inverseTransformOf(relative_center_));
    item->invalidateOpenGLBuffers();
  }
  constraint.setTranslationConstraintType(prev_cons);
}

bool Scene_edit_box_item_priv::applyX(double x, double dirx)
{
  if(x == pool[0])
  {
    if(x+dirx < pool[3])
      return true;
  }
  else
  {
    if(x+dirx > pool[0])
      return true;
  }
  return false;
}

bool Scene_edit_box_item_priv::applyY(double y, double diry)
{
  if(y == pool[1])
  {
    if(y+diry < pool[4])
      return true;
  }
  else
  {
    if(y+diry > pool[1])
      return true;
  }
  return false;
}
bool Scene_edit_box_item_priv::applyZ(double z, double dirz)
{
  if(z == pool[2])
  {
    if(z+dirz < pool[5])
      return true;
  }
  else
  {
    if(z+dirz > pool[2])
      return true;
  }
  return false;
}
