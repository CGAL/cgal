#include "Locally_shortest_path_item.h"
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <sstream>

using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Locally_shortest_path_item_priv Priv;


typedef Scene_surface_mesh_item::Face_graph Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<CGAL::Epick, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

// TODO: update
struct Locally_shortest_path_item::vertex{
  int id;
  double x;
  double y;
  double z;

  CGAL::Point_3<Kernel> position()const
  {
    return  CGAL::Point_3<Kernel>(x,y,z);
  }
  double operator[](int i)
  {
    switch(i)
    {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      return 0;
    }
  }

  template <class P>
  void set(const P& p, int id_)
  {
    x=p.x();
    y=p.y();
    z=p.z();
    id=id_;
  }
};

struct Locally_shortest_path_item_priv{
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

  Locally_shortest_path_item_priv(const CGAL::Three::Scene_interface* scene_interface,
                                  const Scene_surface_mesh_item* sm_item,
                                  Scene_polylines_item* polyline_item,
                                  std::size_t nb_pts,
                                  Locally_shortest_path_item* ebi)
  {
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    ready_to_hl = true;
    scene = scene_interface;
    mesh_item = sm_item;
    spath_item=polyline_item;
    spath_item->polylines.resize(1);
    const Mesh& mesh = *mesh_item->face_graph();
    aabb_tree = Tree(faces(mesh).first,
                     faces(mesh).second,
                     mesh);
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
    vertex_spheres.resize(0);
    normal_spheres.resize(0);
    create_flat_sphere(1.0f, vertex_spheres, normal_spheres,10);


    // TODO: change the default
    if (nb_pts==2)
    {
      //shortest_path
      vertices.resize(2);
      vertices[0].set( mesh.point(*mesh.vertices().begin()), 0 );
      vertices[1].set( mesh.point(*std::next(mesh.vertices().begin())), 1 );
      locations.resize(2);
    }
    else if (nb_pts==4)
    {
      // bezier
      vertices.resize(4);
      vertices[0].set( mesh.point(*mesh.vertices().begin()), 0 );
      vertices[1].set( mesh.point(*std::next(mesh.vertices().begin())), 1 );
      vertices[2].set( mesh.point(*std::next(mesh.vertices().begin(),2)), 2 );
      vertices[3].set( mesh.point(*std::next(mesh.vertices().begin(),3)), 3 );
      locations.resize(4);
    }


    vertex_faces.resize(0);
    normal_faces.resize(0);

    reset_selection();
    last_picked_id = -1;
    last_picked_type = -1;
    QPixmap pix(":/cgal/cursors/resources/rotate_around_cursor.png");
    rotate_cursor = QCursor(pix);

#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
    CGAL::Polygon_mesh_processing::init_geodesic_dual_solver(geodesic_solver, mesh);
#endif
  }

  ~Locally_shortest_path_item_priv(){
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

  bool ready_to_hl;

  CGAL::qglviewer::ManipulatedFrame* frame;
  CGAL::qglviewer::ManipulatedFrame* remodel_frame;
  CGAL::qglviewer::Vec rf_last_pos;
  CGAL::qglviewer::LocalConstraint constraint;
  CGAL::qglviewer::Vec center_;
  CGAL::qglviewer::Vec relative_center_;

  mutable std::vector<Locally_shortest_path_item::vertex> vertices;
  mutable std::vector<PMP::Face_location<Mesh, double>> locations;
  std::vector<Locally_shortest_path_item::vertex*> selected_vertices;

  void reset_selection();
  bool selection_on;
  void picking(int& type, int& id, Viewer_interface *viewer);

  void initializeBuffers(Viewer_interface *viewer)const;

  void computeElements() const;
  void draw_picking(Viewer_interface*);
  void update_points(const QVector3D &dir);

  const Scene_interface* scene;
  const Scene_surface_mesh_item* mesh_item;
  Scene_polylines_item* spath_item;
  Locally_shortest_path_item* item;
  Tree aabb_tree;
  QPoint picked_pixel;
  HL_Primitive hl_type;
  int last_picked_id;
  int last_picked_type;
  QCursor rotate_cursor;
  bool path_invalidated=true;
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
  PMP::Dual_geodesic_solver<double> geodesic_solver;
#endif
};

Locally_shortest_path_item::Locally_shortest_path_item(const CGAL::Three::Scene_interface* scene_interface,
                                                       const Scene_surface_mesh_item *sm_item,
                                                       Scene_polylines_item* polyline_item,
                                                       std::size_t nb_pts)
{
  d = new Locally_shortest_path_item_priv(scene_interface, sm_item, polyline_item, nb_pts, this);

  are_buffers_filled = false;

  for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
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
  contextMenu();
}
QString Locally_shortest_path_item::toolTip() const
{
  std::stringstream ss;
  ss << "Face locations:<br>";
  ss << std::setprecision(17);
  for (auto fl : d->locations)
  {
    ss << "   - " << fl.first << "  (" << fl.second[0] << "," << fl.second[1] << "," << fl.second[2] << ")<br>";
  }
  return QString::fromStdString(ss.str());
}


QMenu* Locally_shortest_path_item::contextMenu()
{
  // disable "Alpha slider" in menu
  QMenu* resMenu = Scene_item::contextMenu();
  bool prop = property("menu_changed").toBool();
  if(!prop)
  {
    setProperty("menu_changed", true);
  }
  return resMenu;
}

void Locally_shortest_path_item::drawSpheres(Viewer_interface *viewer, const QMatrix4x4 f_matrix ) const
{
  GLdouble d_mat[16];
  QMatrix4x4 mv_mat;
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  mv_mat = mv_mat*f_matrix;
  // TODO: must depend on the mesh (and zoom?)
  double radius = 0.01 ;

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

void Locally_shortest_path_item::drawPath() const
{
  if (!d->path_invalidated) return;
  typedef PMP::Edge_location<Mesh, double> Edge_location;
  typedef PMP::Face_location<Mesh, double> Face_location;

  const Mesh& mesh = *d->mesh_item->face_graph();

  if (d->vertices.size()==2)
  {
    std::vector<Edge_location> edge_locations;
    CGAL::Epick::Point_3 src_pt(d->vertices[0].x,d->vertices[0].y,d->vertices[0].z),
                         tgt_pt(d->vertices[1].x,d->vertices[1].y,d->vertices[1].z);

  //TODO store that in the vector vertices
    d->locations[0] = PMP::locate_with_AABB_tree(src_pt, d->aabb_tree, mesh);
    d->locations[1] = PMP::locate_with_AABB_tree(tgt_pt, d->aabb_tree, mesh);

    PMP::locally_shortest_path<double>(d->locations[0], d->locations[1], mesh, edge_locations, d->geodesic_solver);
    d->spath_item->polylines.back().clear();
    d->spath_item->polylines.back().push_back(src_pt);
    for (auto el : edge_locations)
      d->spath_item->polylines.back().push_back(PMP::construct_point(el, mesh));
    d->spath_item->polylines.back().push_back(tgt_pt);
    d->spath_item->setRenderingMode(Wireframe);
    d->spath_item->invalidateOpenGLBuffers();
  }
  else if (d->vertices.size()==4)
  {
    std::vector<Edge_location> edge_locations;
    CGAL::Epick::Point_3 c1_pt(d->vertices[0].x,d->vertices[0].y,d->vertices[0].z),
                         c2_pt(d->vertices[1].x,d->vertices[1].y,d->vertices[1].z),
                         c3_pt(d->vertices[2].x,d->vertices[2].y,d->vertices[2].z),
                         c4_pt(d->vertices[3].x,d->vertices[3].y,d->vertices[3].z);

    //TODO store that in the vector vertices
    d->locations[0] = PMP::locate_with_AABB_tree(c1_pt, d->aabb_tree, mesh);
    d->locations[1] = PMP::locate_with_AABB_tree(c2_pt, d->aabb_tree, mesh);
    d->locations[2] = PMP::locate_with_AABB_tree(c3_pt, d->aabb_tree, mesh);
    d->locations[3] = PMP::locate_with_AABB_tree(c4_pt, d->aabb_tree, mesh);

    PMP::Bezier_segment<Mesh, double> control_points=CGAL::make_array(d->locations[0],
                                                                      d->locations[1],
                                                                      d->locations[2],
                                                                      d->locations[3]);

    std::vector<Face_location> face_locations =
      PMP::recursive_de_Casteljau(mesh, control_points, 8, d->geodesic_solver);

    // TODO: we should connect points with geodesics and not segments
    d->spath_item->polylines.back().clear();
    for (auto fl : face_locations)
      d->spath_item->polylines.back().push_back(PMP::construct_point(fl, mesh));
    d->spath_item->setRenderingMode(Wireframe);
    d->spath_item->invalidateOpenGLBuffers();
  }
  d->path_invalidated=false;
}

void Locally_shortest_path_item::draw(Viewer_interface *viewer) const
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
  drawHl(viewer);

  drawPath();
}

void Locally_shortest_path_item::compute_bbox() const
{
  double xmin=d->vertices[0].x, xmax=xmin;
  double ymin=d->vertices[0].y, ymax=ymin;
  double zmin=d->vertices[0].z, zmax=zmin;
  for(std::size_t i=0; i<d->vertices.size(); ++i)
  {
    xmin=(std::min)(xmin, d->vertices[i].x);
    ymin=(std::min)(ymin, d->vertices[i].y);
    zmin=(std::min)(zmin, d->vertices[i].z);
    xmax=(std::max)(xmax, d->vertices[i].x);
    ymax=(std::max)(ymax, d->vertices[i].y);
    zmax=(std::max)(zmax, d->vertices[i].z);
  }

  setBbox(Scene_item::Bbox(xmin, ymin, zmin,xmax, ymax, zmax));
}

void Locally_shortest_path_item_priv::computeElements() const
{
  vertex_edges.clear();
  vertex_faces.clear();
  normal_faces.clear();
  center_spheres.clear();
  color_edges.clear();
  color_faces.clear();
  color_spheres.clear();

  for (std::size_t i=0; i<vertices.size(); ++i)
  {
    center_spheres.push_back(vertices[i].x);
    center_spheres.push_back(vertices[i].y);
    center_spheres.push_back(vertices[i].z);
  }

  //spheres
  //TODO: check the max value of i that is possible with such a ratio
  for(std::size_t i=0; i<vertices.size(); ++i)
  {
    color_spheres.push_back((20.0*i+10)/255);
    color_spheres.push_back(0);
    color_spheres.push_back(0);
  }
}

Locally_shortest_path_item::~Locally_shortest_path_item()
{
  CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
  viewer->setMouseTracking(false);

  delete d;
}

// Indicate if rendering mode is supported
bool Locally_shortest_path_item::supportsRenderingMode(RenderingMode m) const {
  return m==FlatPlusEdges;
}

Scene_item::ManipulatedFrame*
Locally_shortest_path_item::manipulatedFrame()
{
  return d->frame;
}

void Locally_shortest_path_item::highlight(Viewer_interface *viewer)
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
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].x);
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].y);
      d->hl_vertex.push_back(d->vertices[d->last_picked_id].z);
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
      d->hl_type = Locally_shortest_path_item_priv::VERTEX;
      for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer =
            static_cast<CGAL::Three::Viewer_interface*>(v);
        tc->initializeBuffers(viewer);
      }
      break;
    }
    default:
      d->hl_type = Locally_shortest_path_item_priv::NO_TYPE;
      break;
    }
  }
  else
    clearHL();
  redraw();

  d->ready_to_hl = false;
}

void Locally_shortest_path_item::clearHL()
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
  d->hl_type = Locally_shortest_path_item_priv::NO_TYPE;

  itemChanged();

}
void Locally_shortest_path_item_priv::reset_selection()
{
  selection_on = false;
  CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
  viewer->setManipulatedFrame(frame);
  viewer->setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, CGAL::qglviewer::SELECT);
  constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
  selected_vertices.clear();
}

//intercept events for picking
bool Locally_shortest_path_item::eventFilter(QObject *obj, QEvent *event)
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
      viewer->makeCurrent();
      if(type == 0)
      {
        bool found = false;
        QApplication::setOverrideCursor(Qt::DragMoveCursor);
        CGAL::qglviewer::Vec pos = viewer->camera()->pointUnderPixel(d->picked_pixel, found);
        if(found)
        {
          d->rf_last_pos = pos;
          d->remodel_frame->setPosition(pos);
        }

        d->selection_on = true;
        d->selected_vertices.push_back(&d->vertices[picked]);
        d->constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
        d->remodel_frame->setConstraint(&d->constraint);

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
        CGAL::qglviewer::Vec td(d->remodel_frame->position() - d->rf_last_pos);
        QVector3D dir(td.x, td.y, td.z);
        d->update_points(dir);
        d->rf_last_pos=d->remodel_frame->position();
        d->path_invalidated=true;
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

void Locally_shortest_path_item_priv::draw_picking(Viewer_interface* viewer)
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

  // TODO: radius
  double radius = 0.01 ;
  Tc* tc = item->getTriangleContainer(P_Spheres);
  tc->setFrameMatrix(f_matrix);
  tc->setClipping(false);
  tc->getVao(viewer)->bind();
  tc->getVao(viewer)->program->setAttributeValue("radius", (float)radius);
  tc->getVao(viewer)->release();
  tc->draw(viewer, false);
}

void Locally_shortest_path_item_priv::update_points(const QVector3D &dir)
{
  CGAL::qglviewer::AxisPlaneConstraint::Type prev_cons = constraint.translationConstraintType();
  constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
  for(Locally_shortest_path_item::vertex*  selected_vertex: selected_vertices )
  {
    int id = selected_vertex->id;
    CGAL_assume(id<8 && id >=0);
    double x = selected_vertex->x + dir.x();
    double y = selected_vertex->y + dir.y();
    double z = selected_vertex->z + dir.z();

    // project point onto the mesh
    auto p = aabb_tree.closest_point(CGAL::Epick::Point_3(x,y,z));
    selected_vertex->x = p.x();
    selected_vertex->y = p.y();
    selected_vertex->z = p.z();
  }
  item->invalidateOpenGLBuffers();
  constraint.setTranslationConstraintType(prev_cons);
}

//type : 0 = vertex, 1 = edge, 2 = face
void Locally_shortest_path_item_priv::picking(int& type, int& id, Viewer_interface *viewer)
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

  const auto buffer = read_pixel_as_ubyte_rgba(picked_pixel, viewer, viewer->camera());
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
  viewer->setBackgroundColor(bgColor);
  fbo->release();
  delete fbo;
}

void Locally_shortest_path_item::drawHl(Viewer_interface* viewer)const
{
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

  if(d->hl_type == Locally_shortest_path_item_priv::VERTEX)
  {
    Tc* tc = getTriangleContainer(Priv::S_Spheres);

    tc->setFrameMatrix(f_matrix);
    tc->setMvMatrix(mv_mat);
    tc->setColor(QColor(Qt::yellow));

    // TODO radius
    double radius = 0.01 ;
    tc->setClipping(false);
    tc->getVao(viewer)->bind();
    tc->getVao(viewer)->program->setUniformValue("radius", (float)radius);
    tc->getVao(viewer)->release();
    tc->draw(viewer, true);
  }
}

void Locally_shortest_path_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  setBuffersFilled(false);
  getTriangleContainer(Priv::Spheres)->reset_vbos(ALL);
  getTriangleContainer(Priv::P_Spheres)->reset_vbos(ALL);
}

void Locally_shortest_path_item::computeElements() const
{
  d->computeElements();

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
  setBuffersFilled(true);
}

void Locally_shortest_path_item::initializeBuffers(Viewer_interface *v) const
{
  getTriangleContainer(Priv::Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::P_Spheres)->initializeBuffers(v);
  getTriangleContainer(Priv::P_Spheres)->initializeBuffers(v);

  getTriangleContainer(Priv::Spheres)->setFlatDataSize(d->vertex_spheres.size());
  getTriangleContainer(Priv::Spheres)->setCenterSize(d->center_spheres.size());
  getTriangleContainer(Priv::P_Spheres)->setFlatDataSize(d->vertex_spheres.size());
  getTriangleContainer(Priv::P_Spheres)->setCenterSize(d->center_spheres.size());
}

void Locally_shortest_path_item::connectNewViewer(QObject *o)
{
  Vi* viewer = qobject_cast<Vi*>(o);
  if(!viewer)
    return;
  viewer->setMouseTracking(true);
}
