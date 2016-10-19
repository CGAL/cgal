#include "Scene_edit_box_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>
#include <QGLViewer/manipulatedFrame.h>

struct Scene_edit_box_item::vertex{
  short id;
  CGAL::Point_3<Kernel> position;
  face* face_;
};
struct Scene_edit_box_item::edge{
  short id;
  vertex source;
  vertex target;
  face* face_;
};
struct Scene_edit_box_item::face{
  short id;
  edge edges[4];
  vertex vertices[4];
};

struct Scene_edit_box_item_priv{
  typedef CGAL::Simple_cartesian<double>  Kernel;
  enum VAOs{
    Edges = 0,
    Spheres,
    S_Edges,
    S_Spheres,
    Arrow,
    Faces,
    NumberOfVaos
  };
  enum VBOs{
    VertexEdges = 0,
    VertexSpheres,
    NormalSpheres,
    VertexFaces,
    NormalFaces,
    VertexArrow,
    NormalArrow,
    NumberOfVbos
  };

  Scene_edit_box_item_priv(const CGAL::Three::Scene_interface *scene_interface, Scene_edit_box_item* ebi)
  {
    scene = scene_interface;
    item = ebi;
    CGAL::Three::Scene_item::Bbox bb = scene->bbox();
    double x=(bb.xmin()+bb.xmax())/2;
    double y=(bb.ymin()+bb.ymax())/2;
    double z=(bb.zmin()+bb.zmax())/2;
    center_ = qglviewer::Vec(x,y,z);
    frame = new CGAL::Three::Scene_item::ManipulatedFrame();
    frame->setPosition(center_);
    constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    constraint.setRotationConstraintDirection(qglviewer::Vec(.0,.0,.1));
    frame->setConstraint(&constraint);
    //create the sphere model
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
    for(short i = 0; i< 8; ++i)
    {
      double x,y,z;
      x = ((i/2)%2==0)? scene->bbox().min(0):scene->bbox().max(0);
      y = (i/4==0)? scene->bbox().min(1):scene->bbox().max(1);
      z = (((i+1)/2)%2==1)? scene->bbox().min(2):scene->bbox().max(2);
      vertices[i].position = Kernel::Point_3(x,y,z);
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
    for(short i=0; i<12; ++i)
    {
      edges[i].id = i;
      if(i<4)
      {
        edges[i].source = vertices[i];
        edges[i].target = vertices[i+4];
      }
      else if(i<8)
      {
        edges[i].source = vertices[i];
        edges[i].target = vertices[(i+1)%4 +4];
      }
      else
      {
        edges[i].source = vertices[i%4];
        edges[i].target = vertices[(i+1) %4];
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
    for(short i=0; i<4; ++i)
    {
      faces[0].vertices[i] = vertices[i];
      faces[0].edges[(i+1)%4] = edges[i+7];
    }
    faces[0].id=0;

    for(short i=1; i<4; ++i)
    {
      faces[i].vertices[0] = vertices[i];
      faces[i].vertices[1] = vertices[(i-1)];
      faces[i].vertices[2] = vertices[i+3];
      faces[i].vertices[3] = vertices[i+4];

      faces[i].edges[0] = edges[i+7];
      faces[i].edges[1] = edges[i-1];
      faces[i].edges[2] = edges[i+3];
      faces[i].edges[3] = edges[i];
      faces[i].id = i;
    }

    faces[4].vertices[0] = vertices[0];
    faces[4].vertices[1] = vertices[3];
    faces[4].vertices[2] = vertices[7];
    faces[4].vertices[3] = vertices[4];

    faces[4].edges[0] = edges[0];
    faces[4].edges[1] = edges[3];
    faces[4].edges[2] = edges[7];
    faces[4].edges[3] = edges[4];
    faces[4].id=4;

    for(short i=0; i<4; ++i)
    {
      faces[5].vertices[i] = vertices[i+4];
      faces[5].edges[i] = edges[i+4];
    }
    faces[5].id=5;

    vertex_faces.resize(0);
    normal_faces.resize(0);

    for(short i=0; i<6; ++i)
    {
      for(short j=0; j<4; ++j)
      {
        faces[i].vertices[j].face_ = &faces[i];
        faces[i].edges[j].face_ = &faces[i];
      }
    }
  }

  mutable std::vector<float> vertex_edges;
  mutable std::vector<float> vertex_spheres;
  mutable std::vector<float> normal_spheres;
  mutable std::vector<float> vertex_arrow;
  mutable std::vector<float> normal_arrow;
  mutable std::vector<float> vertex_faces;
  mutable std::vector<float> normal_faces;

  qglviewer::ManipulatedFrame* frame;
  qglviewer::LocalConstraint constraint;
  qglviewer::Vec center_;

  mutable Scene_edit_box_item::vertex vertices[8];
  mutable Scene_edit_box_item::edge edges[12];
  mutable Scene_edit_box_item::face faces[6];

  mutable QOpenGLShaderProgram *program;
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;

  void computeElements() const;

  const CGAL::Three::Scene_interface* scene;
  Scene_edit_box_item* item;
};


Scene_edit_box_item::Scene_edit_box_item(const CGAL::Three::Scene_interface *scene_interface)
  :  Scene_item(NumberOfVbos,NumberOfVaos)

{
  d = new Scene_edit_box_item_priv(scene_interface, this);

  are_buffers_filled = false;
}
QString Scene_edit_box_item::toolTip() const {

  return QString();
}

void Scene_edit_box_item::draw(CGAL::Three::Viewer_interface *viewer) const
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

void Scene_edit_box_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
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
  xmin =d->vertices[0].position.x();
  ymin =d->vertices[0].position.y();
  zmin =d->vertices[0].position.z();
  xmax =d->vertices[6].position.x();
  ymax =d->vertices[6].position.y();
  zmax =d->vertices[6].position.z();
  QVector3D min(xmin, ymin, zmin);
  QVector3D max(xmax, ymax, zmax);
  min = f_matrix*min;
  max = f_matrix*max;
  CGAL::Three::Scene_item::Bbox bb(min.x(),min.y(),min.z(),max.x(),max.y(),max.z());
  _bbox = bb;
}


void
Scene_edit_box_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
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

    item->buffers[VertexEdges].bind();
    program->enableAttributeArray("center");
    program->setAttributeBuffer("center", GL_FLOAT, 0, 3, 3*sizeof(float));
    item->buffers[VertexEdges].release();
    program->disableAttributeArray("radius");
    program->setAttributeValue("radius",1);
    program->disableAttributeArray("color");

    viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
    item->vaos[Spheres]->release();

    program->release();
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

  //edges
  for(short i=0; i<12; ++i)
  {
    if(i<4)
    {
      vertex_edges.push_back(vertices[i].position.x()-center_.x);
      vertex_edges.push_back(vertices[i].position.y()-center_.y);
      vertex_edges.push_back(vertices[i].position.z()-center_.z);

      vertex_edges.push_back(vertices[i+4].position.x()-center_.x);
      vertex_edges.push_back(vertices[i+4].position.y()-center_.y);
      vertex_edges.push_back(vertices[i+4].position.z()-center_.z);
    }
    else if(i<8)
    {
      vertex_edges.push_back(vertices[i].position.x()-center_.x);
      vertex_edges.push_back(vertices[i].position.y()-center_.y);
      vertex_edges.push_back(vertices[i].position.z()-center_.z);

      vertex_edges.push_back(vertices[(i+1)%4 +4].position.x()-center_.x);
      vertex_edges.push_back(vertices[(i+1)%4 +4].position.y()-center_.y);
      vertex_edges.push_back(vertices[(i+1)%4 +4].position.z()-center_.z);
    }
    else
    {
      vertex_edges.push_back(vertices[i%4].position.x()-center_.x);
      vertex_edges.push_back(vertices[i%4].position.y()-center_.y);
      vertex_edges.push_back(vertices[i%4].position.z()-center_.z);

      vertex_edges.push_back(vertices[(i+1) %4].position.x()-center_.x);
      vertex_edges.push_back(vertices[(i+1) %4].position.y()-center_.y);
      vertex_edges.push_back(vertices[(i+1) %4].position.z()-center_.z);
    }
  }
  //faces
  for(short i=0; i<6; ++i)
  {
    push_xyz(vertex_faces, faces[i].vertices[0].position, center_);
    push_xyz(vertex_faces, faces[i].vertices[3].position, center_);
    push_xyz(vertex_faces, faces[i].vertices[2].position, center_);

    push_xyz(vertex_faces, faces[i].vertices[0].position, center_);
    push_xyz(vertex_faces, faces[i].vertices[2].position, center_);
    push_xyz(vertex_faces, faces[i].vertices[1].position, center_);

    for(short j=0; j<6; ++j)
      push_normal(normal_faces, i);
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

CGAL::Three::Scene_item::ManipulatedFrame*
Scene_edit_box_item::manipulatedFrame() { return d->frame; }

double Scene_edit_box_item::point(short i, short j) const
{
  QVector3D pos(d->vertices[i].position.x()-d->center_.x,
                d->vertices[i].position.y()-d->center_.y,
                d->vertices[i].position.z()-d->center_.z);
  QMatrix4x4 f_matrix;
  for (int k=0; k<16; ++k){
    f_matrix.data()[k] = (float)d->frame->matrix()[k];
  }
  return (f_matrix*pos)[j];
}
