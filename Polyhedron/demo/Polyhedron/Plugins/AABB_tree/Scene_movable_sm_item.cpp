#include <QApplication>
#include "Scene_movable_sm_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
typedef CGAL::Three::Triangle_container Tri;
typedef CGAL::Three::Viewer_interface VI;
struct Scene_movable_sm_item_priv
{
  Scene_movable_sm_item_priv(const CGAL::qglviewer::Vec& pos,SMesh* sm,
                             const QString name, Scene_movable_sm_item *parent)
    : frame(new CGAL::Three::Scene_item::ManipulatedFrame()),
      facegraph(sm),
      center_(pos),
      item_name(name)
  {
    item = parent;
    const CGAL::qglviewer::Vec offset = static_cast<Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    frame->setPosition(pos+offset);
    item->setTriangleContainer(0, new Triangle_container(VI::PROGRAM_WITH_LIGHT,
                                                         false));
    item->setEdgeContainer(0, new Edge_container(VI::PROGRAM_NO_SELECTION,
                                                 false));
  }
  ~Scene_movable_sm_item_priv()
  {
    delete frame;
  }
  void initialize_buffers(Viewer_interface *viewer) const;
  void compute_elements() const;
  enum VAOs {
    Edges=0,
    NbOfVaos
  };
  enum VBOs {
    Vertices = 0,
    NbOfVbos
  };

  CGAL::qglviewer::ManipulatedFrame* frame;
  SMesh* facegraph;
  CGAL::qglviewer::Vec center_;
  Scene_movable_sm_item *item;
  QMatrix4x4 f_matrix;
  const QString item_name;

  mutable QOpenGLShaderProgram *program;
  mutable std::vector<float> flat_normals;
  mutable std::vector<float> flat_vertices;
  mutable std::vector<float> edges_vertices;

};

Scene_movable_sm_item::Scene_movable_sm_item(const CGAL::qglviewer::Vec& pos, SMesh* sm,
                                             const QString name)
{
  d = new Scene_movable_sm_item_priv(pos,sm, name, this);
}


void Scene_movable_sm_item_priv::initialize_buffers(CGAL::Three::Viewer_interface *viewer = nullptr) const
{
  item->getTriangleContainer(0)->initializeBuffers(viewer);
  item->getTriangleContainer(0)->setFlatDataSize(flat_vertices.size());
  item->getEdgeContainer(0)->initializeBuffers(viewer);
  item->getEdgeContainer(0)->setFlatDataSize(edges_vertices.size());
  flat_vertices.resize(0);
  flat_normals .resize(0);
  edges_vertices.resize(0);
  flat_vertices.shrink_to_fit();
  flat_normals.shrink_to_fit();
  edges_vertices.shrink_to_fit();

  item->are_buffers_filled = true;
}


void Scene_movable_sm_item_priv::compute_elements() const
{
  typedef EPICK::Point_3 Point;
  QApplication::setOverrideCursor(Qt::WaitCursor);

  SMesh::Property_map<face_descriptor, EPICK::Vector_3 > fnormals =
      facegraph->add_property_map<face_descriptor, EPICK::Vector_3 >("f:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*facegraph,fnormals);
  const CGAL::qglviewer::Vec o = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  EPICK::Vector_3 offset(o.x, o.y, o.z);
  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
      facegraph->points();
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;
  typedef CGAL::Buffer_for_vao<float, unsigned int> CPF;
  flat_vertices.clear();
  flat_normals.clear();
  edges_vertices.clear();
  //faces
  for(face_descriptor fd : faces(*facegraph))
  {
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *facegraph),*facegraph))
    {
      Point p = positions[source(hd, *facegraph)] + offset;
      EPICK::Point_3 pc(p.x() - center_.x,
                        p.y() - center_.y,
                        p.z() - center_.z);
      CPF::add_point_in_buffer(pc, flat_vertices);
      EPICK::Vector_3 n = fnormals[fd];
      CPF::add_normal_in_buffer(n, flat_normals);
    }
  }
  //edges
  for(edge_descriptor ed : edges(*facegraph))
  {
    Point p = positions[source(ed, *facegraph)] + offset;
    EPICK::Point_3 pc(p.x() - center_.x,
                      p.y() - center_.y,
                      p.z() - center_.z);
    CPF::add_point_in_buffer(pc, edges_vertices);
    p = positions[target(ed, *facegraph)] + offset;
    pc=EPICK::Point_3(p.x() - center_.x,
                      p.y() - center_.y,
                      p.z() - center_.z);
    CPF::add_point_in_buffer(pc, edges_vertices);
  }



  item->getTriangleContainer(0)->allocate(Tri::Flat_vertices, flat_vertices.data(),
                                          static_cast<int>(flat_vertices.size()*sizeof(float)));
  item->getTriangleContainer(0)->allocate(Tri::Flat_normals, flat_normals.data(),
                                          static_cast<int>(flat_normals.size()*sizeof(float)));
  item->getEdgeContainer(0)->allocate(Tri::Flat_vertices, edges_vertices.data(),
                                      static_cast<int>(edges_vertices.size()*sizeof(float)));
  QApplication::restoreOverrideCursor();
}

void Scene_movable_sm_item::computeElements() const
{
  d->compute_elements();
  compute_bbox();
}
void Scene_movable_sm_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer) && viewer->context()->isValid())
    initGL(viewer);
  if(!are_buffers_filled)
    d->initialize_buffers(viewer);
  getTriangleContainer(0)->setColor(color());
  getTriangleContainer(0)->setAlpha(alpha());
  getTriangleContainer(0)->setFrameMatrix(d->f_matrix);

  getTriangleContainer(0)->draw(viewer, true);
}

void Scene_movable_sm_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer) && viewer->context()->isValid())
    initGL(viewer);
  if(!are_buffers_filled)
    d->initialize_buffers(viewer);
  getEdgeContainer(0)->setColor(color());
  getEdgeContainer(0)->setFrameMatrix(d->f_matrix);

  getEdgeContainer(0)->draw(viewer, true);
}

QString Scene_movable_sm_item::toolTip() const {
  return QObject::tr("<p>Manipulatable representation of <b>%1</b></p>"
                     "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />")
      .arg(d->item_name);
}

void
Scene_movable_sm_item::compute_bbox() const {
  SMesh::Property_map<vertex_descriptor, Point_3> pprop = d->facegraph->points();
  CGAL::Bbox_3 bbox ;

  for(vertex_descriptor vd :vertices(*d->facegraph))
  {
    bbox = bbox + pprop[vd].bbox();
  }
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax());
  is_bbox_computed = true;
}

Scene_item::Bbox Scene_movable_sm_item::bbox() const {
  if(!is_bbox_computed)
    compute_bbox();
  is_bbox_computed = true;
  return _bbox;
}


void Scene_movable_sm_item::invalidateOpenGLBuffers()
{
  d->compute_elements();
  is_bbox_computed = false;
  are_buffers_filled = false;
}
CGAL::Three::Scene_item::ManipulatedFrame* Scene_movable_sm_item::manipulatedFrame() { return d->frame; }
const CGAL::qglviewer::Vec& Scene_movable_sm_item::center() const { return d->center_; }
Scene_movable_sm_item::~Scene_movable_sm_item() { delete d; Q_EMIT killed(); }
void Scene_movable_sm_item::setFMatrix(const GLdouble matrix[])
{
  for (int i=0; i<16; ++i)
    d->f_matrix.data()[i] = (float)matrix[i];
}

SMesh *Scene_movable_sm_item::getFaceGraph()
{
  return d->facegraph;
}
