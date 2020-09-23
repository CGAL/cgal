#include "Scene_edit_box_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/constraint.h>
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLShaderProgram>
#include <QCursor>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Scene_edit_box_item_priv Priv;

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
  enum Face_containers{
    Faces = 0,
    S_Faces,
    Spheres,
    S_Spheres,
    P_Spheres,
    P_Faces,
    Nbf
  };

  enum Line_containers{
    Edges = 0,
    S_Edges,
    P_Edges,
    Nbe
  };

  enum HL_Primitive{
    VERTEX=0,
    EDGE,
    FACE,
    NO_TYPE
  };

  Scene_edit_box_item_priv(const Scene_interface *scene_interface, Scene_edit_box_item* ebi)
  {
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    ready_to_hl = true;
    scene = scene_interface;
    item = ebi;
    selection_on = false;
    Scene_item::Bbox bb = scene->bbox();
    double x=(bb.xmin()+bb.xmax())/2;
    double y=(bb.ymin()+bb.ymax())/2;
    double z=(bb.zmin()+bb.zmax())/2;
    center_ = CGAL::qglviewer::Vec(x,y,z);
    relative_center_ = CGAL::qglviewer::Vec(0,0,0);
    remodel_frame = new Scene_item::ManipulatedFrame();
    remodel_frame->setTranslationSensitivity(1.0);
    frame = new Scene_item::ManipulatedFrame();
    frame->setPosition(center_+offset);
    frame->setSpinningSensitivity(100.0); //forbid spinning
    constraint.setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::AXIS);
    constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
    constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(.0,.0,.1));
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


  CGAL::qglviewer::ManipulatedFrame* frame;
  CGAL::qglviewer::ManipulatedFrame* remodel_frame;
  CGAL::qglviewer::Vec rf_last_pos;
  CGAL::qglviewer::LocalConstraint constraint;
  CGAL::qglviewer::Vec center_;
  CGAL::qglviewer::Vec relative_center_;

  mutable Scene_edit_box_item::vertex vertices[8];
  mutable Scene_edit_box_item::edge edges[12];
  mutable Scene_edit_box_item::face faces[6];

  std::vector<Scene_edit_box_item::vertex*> selected_vertices;
  void reset_selection();
  bool selection_on;
  void picking(int& type, int& id, Viewer_interface *viewer);

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

Scene_edit_box_item::Scene_edit_box_item()
{
  d = NULL;
}
Scene_edit_box_item::Scene_edit_box_item(const Scene_interface *scene_interface)
{
  d = new Scene_edit_box_item_priv(scene_interface, this);

  are_buffers_filled = false;

  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    v->setMouseTracking(true);
  }
  connect(Three::mainWindow(), SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));

  setTriangleContainer(Priv::P_Faces  , new Tc(Vi::PROGRAM_NO_SELECTION, false));
  setTriangleContainer(Priv::Faces , new Tc(Vi::PROGRAM_WITH_LIGHT, false));
  setTriangleContainer(Priv::S_Faces , new Tc(Vi::PROGRAM_WITH_LIGHT, false));
  setTriangleContainer(Priv::Spheres  , new Tc(Vi::PROGRAM_SPHERES, false));
  setTriangleContainer(Priv::S_Spheres, new Tc(Vi::PROGRAM_SPHERES, false));
  setTriangleContainer(Priv::P_Spheres, new Tc(Vi::PROGRAM_DARK_SPHERES, false));


  for(int i=Priv::Nbe-1; i>=0; --i)
  {
    setEdgeContainer(i,
                     new Ec(Three::mainViewer()->isOpenGL_4_3()
                            ? Vi::PROGRAM_SOLID_WIREFRAME
                            : Vi::PROGRAM_NO_SELECTION,
                            false));
  }
}
QString Scene_edit_box_item::toolTip() const {

  return QString();
}

void Scene_edit_box_item::drawSpheres(Viewer_interface *viewer, const QMatrix4x4 f_matrix ) const
{
  GLdouble d_mat[16];
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  double radius =std::sqrt(
      (point(6,0) - point(0,0)) * (point(6,0) - point(0,0)) +
      (point(6,1) - point(0,1)) * (point(6,1) - point(0,1)) +
      (point(6,2) - point(0,2)) * (point(6,2) - point(0,2))) *0.02 ;

  Tc* tc = getTriangleContainer(Priv::Spheres);
  tc->setFrameMatrix(f_matrix);
  tc->setMvMatrix(mv_mat);
  tc->setClipping(false);
  tc->getVao(viewer)->bind();
  tc->getVao(viewer)->program->setAttributeValue("radius",radius);
  tc->getVao(viewer)->release();
  tc->setColor(QColor(Qt::red));
  tc->draw(viewer, true);
}

void Scene_edit_box_item::draw(Viewer_interface *viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }

  drawSpheres(viewer, f_matrix);

}

void Scene_edit_box_item::drawEdges(Viewer_interface* viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }
  Ec* ec = getEdgeContainer(Priv::Edges);
  if(viewer->isOpenGL_4_3())
  {
    QVector2D vp(viewer->width(), viewer->height());
    ec->setViewport(vp);
    ec->setWidth(6.0f);
  }
  ec->setClipping(false);
  ec->setFrameMatrix(f_matrix);
  ec->setColor(QColor(Qt::black));
  ec->draw(viewer, true);

  if(renderingMode() == Wireframe)
  {
    drawSpheres(viewer, f_matrix);
  }
  drawHl(viewer);
  drawTransparent(viewer);
}

void Scene_edit_box_item::compute_bbox() const
{
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();


  QVector3D vmin(d->pool[0], d->pool[1], d->pool[2]);
  QVector3D vmax(d->pool[3], d->pool[4], d->pool[5]);

  for(int i=0; i< 3; ++i)
  {
    vmin[i] += d->frame->translation()[i]-d->center_[i]-offset[i];
    vmax[i] += d->frame->translation()[i]-d->center_[i]-offset[i];
  }

  setBbox(Scene_item::Bbox(vmin.x(),vmin.y(),vmin.z(),vmax.x(),vmax.y(),vmax.z()));
}




void push_xyz(std::vector<float> &v,
              const Scene_edit_box_item::Kernel::Point_3& p,
              CGAL::qglviewer::Vec center_ = CGAL::qglviewer::Vec(0,0,0))
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
  CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
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
  CGAL::qglviewer::Vec pos(d->vertices[i].position().x()-d->center_.x,
                     d->vertices[i].position().y()-d->center_.y,
                     d->vertices[i].position().z()-d->center_.z);
  return (d->frame->inverseCoordinatesOf(pos))[j];
}

void Scene_edit_box_item::highlight(Viewer_interface *viewer)
{
  d->ready_to_hl = true;
  viewer->makeCurrent();
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
      Tc* tc = getTriangleContainer(Priv::S_Spheres);
      tc->reset_vbos(ALL);
      tc->allocate(
            Tc::Flat_vertices,
            d->vertex_spheres.data(),
            static_cast<int>(d->vertex_spheres.size()*sizeof(float)));

      tc->allocate(
            Tc::Flat_normals,
            d->normal_spheres.data(),
            static_cast<int>(d->normal_spheres.size()*sizeof(float)));

      tc->allocate(
            Tc::Facet_centers,
            d->hl_vertex.data(),
            static_cast<int>(d->hl_vertex.size()*sizeof(float)));

      tc->setFlatDataSize(d->vertex_spheres.size());
      tc->setCenterSize(d->hl_vertex.size());
      //draw
      d->hl_type = Scene_edit_box_item_priv::VERTEX;
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer =
            static_cast<CGAL::Three::Viewer_interface*>(v);
        tc->initializeBuffers(viewer);
      }
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
      Ec* ec = getEdgeContainer(Priv::S_Edges);
      ec->reset_vbos(ALL);
      ec->allocate(
            Ec::Vertices,
            d->hl_vertex.data(),
            static_cast<GLsizei>(d->hl_vertex.size()*sizeof(float)));
      ec->setFlatDataSize(d->hl_vertex.size());
      //draw
      d->hl_type = Scene_edit_box_item_priv::EDGE;
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer =
            static_cast<CGAL::Three::Viewer_interface*>(v);
        ec->initializeBuffers(viewer);
      }
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
      Tc* tc = getTriangleContainer(Priv::S_Faces);
      tc->reset_vbos(ALL);
      tc->allocate(
            Tc::Flat_vertices,
            d->hl_vertex.data(),
            static_cast<int>(d->hl_vertex.size()*sizeof(float)));

      tc->allocate(
            Tc::Flat_normals,
            d->hl_normal.data(),
            static_cast<int>(d->hl_normal.size()*sizeof(float)));
      tc->setFlatDataSize(d->hl_vertex.size());
      //draw
      d->hl_type = Scene_edit_box_item_priv::FACE;
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer =
            static_cast<CGAL::Three::Viewer_interface*>(v);
        tc->initializeBuffers(viewer);
      }
      break;
    }
    default:
      d->hl_type = Scene_edit_box_item_priv::NO_TYPE;
      break;
    }
  }
  else
    clearHL();
  redraw();

  d->ready_to_hl = false;
}

void Scene_edit_box_item::clearHL()
{
  Viewer_interface* viewer = dynamic_cast<Viewer_interface*>(*CGAL::QGLViewer::QGLViewerPool().begin());
  viewer->makeCurrent();
  d->hl_normal.clear();
  d->hl_vertex.clear();

  Tc* tc = getTriangleContainer(Priv::S_Spheres);
  tc->reset_vbos(ALL);
  tc->allocate(Tc::Flat_vertices, d->vertex_spheres.data(),
               static_cast<int>(d->vertex_spheres.size()*sizeof(float)));
  tc->allocate(Tc::Flat_normals,
               d->normal_spheres.data(),
               static_cast<int>(d->normal_spheres.size()*sizeof(float)));

  tc->allocate(Tc::Facet_centers, nullptr, 0);

  tc->initializeBuffers(viewer);
  tc->setFlatDataSize(0);
  tc->setCenterSize(0);
  //draw
  Ec* ec = getEdgeContainer(Priv::S_Edges);
  ec->reset_vbos(ALL);
  ec->allocate(Ec::Vertices, nullptr, 0);
  ec->initializeBuffers(viewer);
  ec->setFlatDataSize(0);

  tc = getTriangleContainer(Priv::S_Faces);
  tc->reset_vbos(ALL);
  tc->allocate(Tc::Flat_vertices, nullptr, 0);
  tc->allocate(Tc::Flat_normals, nullptr, 0);
  tc->initializeBuffers(viewer);
  tc->setFlatDataSize(0);
  d->hl_type = Scene_edit_box_item_priv::NO_TYPE;

  itemChanged();

}
void Scene_edit_box_item_priv::reset_selection()
{
  selection_on = false;
  CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
  viewer->setManipulatedFrame(frame);
  viewer->setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, CGAL::qglviewer::SELECT);
  constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
  selected_vertices.clear();
}

//intercept events for picking
bool Scene_edit_box_item::eventFilter(QObject *obj, QEvent *event)
{
  if(!visible())
    return false;
  Viewer_interface* viewer = qobject_cast<Vi*>(obj);
  if(!viewer)
    return false;
  if(event->type() == QEvent::MouseButtonPress)
  {
    QMouseEvent* e = static_cast<QMouseEvent*>(event);
    if(e->modifiers() == Qt::NoModifier)
    {
      //pick
      int type, picked;
      d->picked_pixel = e->pos();
      d->picking(type, picked, viewer);
      if(type !=-1)
      {
        bool found = false;
        QApplication::setOverrideCursor(Qt::DragMoveCursor);
        CGAL::qglviewer::Vec pos = viewer->camera()->pointUnderPixel(d->picked_pixel, found);
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
          d->constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
          d->remodel_frame->setConstraint(&d->constraint);
        }
        else if(type == 1)
        {
          d->selected_vertices.push_back(d->edges[picked].source);
          d->selected_vertices.push_back(d->edges[picked].target);
          Kernel::Point_3 s(d->edges[picked].source->position()), t(d->edges[picked].target->position());

          CGAL::qglviewer::Vec normal(t.x()-s.x(), t.y()-s.y(), t.z()-s.z());
          d->constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::PLANE);
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
          d->constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::AXIS);
          d->constraint.setTranslationConstraintDirection(CGAL::qglviewer::Vec(n.x(), n.y(), n.z()));
        }

        viewer->setManipulatedFrame(d->remodel_frame);
        viewer->setMouseBinding(
                    Qt::NoModifier,
                    Qt::LeftButton,
                    CGAL::qglviewer::FRAME,
                    CGAL::qglviewer::TRANSLATE);
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
    if(e->modifiers() == Qt::NoModifier)
    {
      if(d->selection_on)
      {
        d->remodel_frame->setOrientation(d->frame->orientation());
        CGAL::qglviewer::Vec td(d->remodel_frame->transformOf(d->remodel_frame->position()
                                                              - d->rf_last_pos));
        QVector3D dir(td.x, td.y, td.z);
        d->remodel_box(dir);
      }
      d->ready_to_hl= true;
      d->picked_pixel = e->pos();
      QTimer::singleShot(0, this,
                         [this, viewer](){
        highlight(viewer);
      });
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
    viewer->setMouseBinding(
                Qt::NoModifier,
                Qt::LeftButton,
                CGAL::qglviewer::CAMERA,
                CGAL::qglviewer::ROTATE);
  }
  else if(event->type() == QEvent::KeyRelease)
  {
     QKeyEvent* e = static_cast<QKeyEvent*>(event);
     if(e->key() == Qt::Key_Control)
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
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;


  if(item->renderingMode() == FlatPlusEdges)
  {
    Tc* tc = item->getTriangleContainer(P_Faces);
    tc->setFrameMatrix(f_matrix);
    tc->setClipping(false);
    tc->draw(viewer, false);
  }
  double radius =std::sqrt(
      (item->point(6,0) - item->point(0,0)) * (item->point(6,0) - item->point(0,0)) +
      (item->point(6,1) - item->point(0,1)) * (item->point(6,1) - item->point(0,1)) +
      (item->point(6,2) - item->point(0,2)) * (item->point(6,2) - item->point(0,2))) *0.02 ;
  Tc* tc = item->getTriangleContainer(P_Spheres);
  tc->setFrameMatrix(f_matrix);
  tc->setClipping(false);
  tc->getVao(viewer)->bind();
  tc->getVao(viewer)->program->setAttributeValue("radius", (float)radius);
  tc->getVao(viewer)->release();
  tc->draw(viewer, false);

  Ec* ec = item->getEdgeContainer(P_Edges);
  if(viewer->isOpenGL_4_3())
  {

    QVector2D vp(viewer->width(), viewer->height());
    ec->setViewport(vp);
    ec->setWidth(6.0f);
  }
  ec->setFrameMatrix(f_matrix);
  ec->setClipping(false);
  ec->draw(viewer, false);
}

void Scene_edit_box_item_priv::remodel_box(const QVector3D &dir)
{
  CGAL::qglviewer::AxisPlaneConstraint::Type prev_cons = constraint.translationConstraintType();
  constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
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
  viewer->makeCurrent();
  type = -1;
  id = -1;
  int deviceWidth = viewer->camera()->screenWidth();
  int deviceHeight = viewer->camera()->screenHeight();
  QOpenGLFramebufferObject* fbo = new QOpenGLFramebufferObject(deviceWidth, deviceHeight,QOpenGLFramebufferObject::Depth);
  fbo->bind();
  viewer->glEnable(GL_DEPTH_TEST);
  viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  QColor bgColor(viewer->backgroundColor());
  //draws the image in the fbo
  viewer->setBackgroundColor(::Qt::white);
  draw_picking(viewer);

  int rowLength = deviceWidth * 4; // data asked in RGBA,so 4 bytes.
  const static int dataLength = rowLength * deviceHeight;
  GLubyte* buffer = new GLubyte[dataLength];
  // Qt uses upper corner for its origin while GL uses the lower corner.
  viewer->glReadPixels(picked_pixel.x(), deviceHeight-1-picked_pixel.y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
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
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;

  if(d->hl_type == Scene_edit_box_item_priv::VERTEX)
  {
    Tc* tc = getTriangleContainer(Priv::S_Spheres);

    tc->setFrameMatrix(f_matrix);
    tc->setMvMatrix(mv_mat);
    tc->setColor(QColor(Qt::yellow));

    double radius =std::sqrt(
        (point(6,0) - point(0,0)) * (point(6,0) - point(0,0)) +
        (point(6,1) - point(0,1)) * (point(6,1) - point(0,1)) +
        (point(6,2) - point(0,2)) * (point(6,2) - point(0,2))) *0.02 ;
    tc->setClipping(false);
    tc->getVao(viewer)->bind();
    tc->getVao(viewer)->program->setUniformValue("radius", (float)radius);
    tc->getVao(viewer)->release();
    tc->draw(viewer, true);
  }
  else if(d->hl_type == Scene_edit_box_item_priv::EDGE)
  {
    Ec* ec = getEdgeContainer(Priv::S_Edges);
    if(viewer->isOpenGL_4_3())
    {
      QVector2D vp(viewer->width(), viewer->height());
      ec->setViewport(vp);
      ec->setWidth(6.0f);
    }
    ec->setClipping(false);
    ec->setFrameMatrix(f_matrix);
    ec->setColor(QColor(Qt::yellow));
    ec->draw(viewer, true);

  }
  else if(d->hl_type == Scene_edit_box_item_priv::FACE)
  {
    viewer->glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
    viewer->glEnable(GL_BLEND);
    viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Tc* tc = getTriangleContainer(Priv::S_Faces);
    tc->setMvMatrix(mv_mat);
    tc->setFrameMatrix(f_matrix);
    tc->setClipping(false);

    tc->setColor(QColor(Qt::yellow));
    tc->setAlpha(0.5);
    tc->draw(viewer, true);
    viewer->glPolygonOffset(offset_factor, offset_units);
    viewer->glDisable(GL_BLEND);
  }
}
void Scene_edit_box_item::drawTransparent(CGAL::Three::Viewer_interface*viewer)const
{
  if(renderingMode() != FlatPlusEdges)
    return;
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  QMatrix4x4 f_matrix;
  for (int i=0; i<16; ++i){
    f_matrix.data()[i] = (float)d->frame->matrix()[i];
  }

  GLdouble d_mat[16];
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;

  viewer->glEnable(GL_BLEND);
  viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  Tc* tc = getTriangleContainer(Priv::Faces);
  tc->setMvMatrix(mv_mat);
  tc->setFrameMatrix(f_matrix);
  tc->setClipping(false);
  tc->setColor(QColor(128,128,128,128));
  tc->setAlpha(0.5);
  tc->draw(viewer, true);
  viewer->glDisable(GL_BLEND);
}

void Scene_edit_box_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  setBuffersFilled(false);
  getTriangleContainer(Priv::Faces)->reset_vbos(ALL);
  getTriangleContainer(Priv::P_Faces)->reset_vbos(ALL);
  getTriangleContainer(Priv::Spheres)->reset_vbos(ALL);
  getTriangleContainer(Priv::P_Spheres)->reset_vbos(ALL);
  getEdgeContainer(Priv::Edges)->reset_vbos(ALL);
  getEdgeContainer(Priv::P_Edges)->reset_vbos(ALL);
}

void Scene_edit_box_item::computeElements() const
{
  d->computeElements();
  getEdgeContainer(Priv::Edges)->allocate(
        Ec::Vertices,
        d->vertex_edges.data(),
        static_cast<GLsizei>(d->vertex_edges.size()*sizeof(float)));
  Ec* ec = getEdgeContainer(Priv::P_Edges);
  ec->allocate(
        Ec::Vertices,
        d->vertex_edges.data(),
        static_cast<GLsizei>(d->vertex_edges.size()*sizeof(float)));

  ec->allocate(
        Ec::Colors,
        d->color_edges.data(),
        static_cast<GLsizei>(d->color_edges.size()*sizeof(float)));


  Tc* tc = getTriangleContainer(Priv::Spheres);
  tc->allocate(
        Tc::Flat_vertices,
        d->vertex_spheres.data(),
        static_cast<int>(d->vertex_spheres.size()*sizeof(float)));

  tc->allocate(
        Tc::Flat_normals,
        d->normal_spheres.data(),
        static_cast<int>(d->normal_spheres.size()*sizeof(float)));
  tc->allocate(
        Tc::Facet_centers,
        d->center_spheres.data(),
        static_cast<int>(d->center_spheres.size()*sizeof(float)));

  tc = getTriangleContainer(Priv::P_Spheres);
  tc->allocate(
        Tc::Flat_vertices,
        d->vertex_spheres.data(),
        static_cast<int>(d->vertex_spheres.size()*sizeof(float)));

  tc->allocate(
        Tc::Facet_centers,
        d->center_spheres.data(),
        static_cast<int>(d->center_spheres.size()*sizeof(float)));

  tc->allocate(
        Tc::FColors,
        d->color_spheres.data(),
        static_cast<int>(d->color_spheres.size()*sizeof(float)));

  tc = getTriangleContainer(Priv::Faces);
  tc->allocate(
        Tc::Flat_vertices,
        d->vertex_faces.data(),
        static_cast<int>(d->vertex_faces.size()*sizeof(float)));

  tc->allocate(
        Tc::Flat_normals,
        d->normal_faces.data(),
        static_cast<int>(d->normal_faces.size()*sizeof(float)));
  tc = getTriangleContainer(Priv::P_Faces);
  tc->allocate(
        Tc::Flat_vertices,
        d->vertex_faces.data(),
        static_cast<int>(d->vertex_faces.size()*sizeof(float)));
  tc->allocate(
        Tc::FColors,
        d->color_faces.data(),
        static_cast<int>(d->color_faces.size()*sizeof(float)));
  setBuffersFilled(true);
}

void Scene_edit_box_item::initializeBuffers(Viewer_interface *v) const
{

  getTriangleContainer(Priv::Faces)->initializeBuffers(v);
  getTriangleContainer(Priv::P_Faces)->initializeBuffers(v);

  getTriangleContainer(Priv::Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::P_Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::P_Spheres)->initializeBuffers(v);

  getEdgeContainer(Priv::Edges)->initializeBuffers(v);
  getEdgeContainer(Priv::P_Edges)->initializeBuffers(v);


  getTriangleContainer(Priv::Faces)->setFlatDataSize(d->vertex_faces.size());
  getTriangleContainer(Priv::P_Faces)->setFlatDataSize(d->vertex_faces.size());

  getTriangleContainer(Priv::Spheres)->setFlatDataSize(d->vertex_spheres.size());
  getTriangleContainer(Priv::Spheres)->setCenterSize(d->center_spheres.size());
  getTriangleContainer(Priv::P_Spheres)->setFlatDataSize(d->vertex_spheres.size());
  getTriangleContainer(Priv::P_Spheres)->setCenterSize(d->center_spheres.size());

  getEdgeContainer(Priv::Edges)->setFlatDataSize(d->vertex_edges.size());
  getEdgeContainer(Priv::P_Edges)->setFlatDataSize(d->vertex_edges.size());
}

void Scene_edit_box_item::connectNewViewer(QObject *o)
{
  Vi* viewer = qobject_cast<Vi*>(o);
  if(!viewer)
    return;
  viewer->setMouseTracking(true);
}
