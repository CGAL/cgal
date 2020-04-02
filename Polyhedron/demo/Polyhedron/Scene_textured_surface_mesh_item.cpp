#include "Scene_textured_surface_mesh_item.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>
#include <QApplication>
#include <QObject>

typedef EPICK  ::Point_3 Point;
using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;

struct Scene_textured_surface_mesh_item_priv
{
  Scene_textured_surface_mesh_item_priv(Scene_textured_surface_mesh_item* parent)
    :sm(new SMesh)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    umap = sm->add_property_map<halfedge_descriptor, float>("h:u", 0.0f).first;
    vmap = sm->add_property_map<halfedge_descriptor, float>("h:v", 0.0f).first;
  }
  Scene_textured_surface_mesh_item_priv(const SMesh& p, Scene_textured_surface_mesh_item* parent)
    : sm(new SMesh(p))

  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    umap = sm->add_property_map<halfedge_descriptor, float>("h:u", 0.0f).first;
    vmap = sm->add_property_map<halfedge_descriptor, float>("h:v", 0.0f).first;
  }
  Scene_textured_surface_mesh_item_priv(SMesh* const p,Scene_textured_surface_mesh_item* parent)
    :sm(p)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    umap = sm->add_property_map<halfedge_descriptor, float>("h:u", 0.0f).first;
    vmap = sm->add_property_map<halfedge_descriptor, float>("h:v", 0.0f).first;
  }

  ~Scene_textured_surface_mesh_item_priv()
  {
    delete sm;
  }

  void compute_normals_and_vertices(void) const;

  SMesh* sm;
  ::Texture texture;
  SMesh::Property_map<halfedge_descriptor, float> umap;
  SMesh::Property_map<halfedge_descriptor, float> vmap;

  //[Px][Py][Pz][Nx][Ny][Nz][u][v]
  mutable std::vector<float> faces_buffer;
  //[Px][Py][Pz][u][v]
  mutable std::vector<float> edges_buffer;
  //[Px][Py][Pz]
  mutable std::vector<float> border_edges_buffer;
  mutable std::size_t nb_face_verts;
  mutable std::size_t nb_edge_verts;
  mutable std::size_t nb_border_verts;
  bool are_buffers_filled;
  bool are_border_filled;

  Scene_textured_surface_mesh_item* item;
  typedef Scene_textured_surface_mesh_item I;
};
void
Scene_textured_surface_mesh_item_priv::compute_normals_and_vertices(void) const
{
  faces_buffer.resize(0);

  typedef boost::graph_traits<SMesh>::face_iterator face_iterator;
  typedef boost::graph_traits<SMesh>::face_iterator face_iterator;
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();

  //Faces
  SMesh::Property_map<vertex_descriptor, Point> positions =
      sm->points();
  SMesh::Property_map<face_descriptor, EPICK::Vector_3 > fnormals =
      sm->add_property_map<face_descriptor, EPICK::Vector_3 >("f:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*sm,fnormals);

  for(face_iterator f = faces(*sm).begin();
      f != faces(*sm).end();
      ++f)
  {

    SMesh::Halfedge_around_face_circulator he(halfedge(*f, *sm), *sm);
    SMesh::Halfedge_around_face_circulator end = he;
    CGAL_For_all(he,end)
    {
      //position [3]
      const Point& p = get(positions, target(*he, *sm));
      faces_buffer.push_back(p.x() + offset.x);
      faces_buffer.push_back(p.y() + offset.y);
      faces_buffer.push_back(p.z() + offset.z);
      //normals [3]
      const EPICK::Vector_3& n = get(fnormals, face(*he, *sm));
      faces_buffer.push_back(n[0]);
      faces_buffer.push_back(n[1]);
      faces_buffer.push_back(n[2]);
      //uvs [2]
      const float u = get(umap, *he);
      const float v = get(vmap, *he);
      faces_buffer.push_back(u);
      faces_buffer.push_back(v);
    }
  }

  //Edges
  typedef EPICK::Point_3                Point;
  typedef SMesh::Edge_iterator        Edge_iterator;

  Edge_iterator he;

  for(he = edges(*sm).begin();
      he != edges(*sm).end();
      ++he)
  {

    //position [3]
      const Point& a = sm->point(target(*he, *sm));
      const Point& b = sm->point(source(*he, *sm));
      edges_buffer.push_back(a.x() + offset.x);
      edges_buffer.push_back(a.y() + offset.y);
      edges_buffer.push_back(a.z() + offset.z);
    //uvs [2]
      float u = get(umap, halfedge(*he, *sm));
      float v = get(vmap, halfedge(*he, *sm));

      edges_buffer.push_back(u);
      edges_buffer.push_back(v);
      //position [3]
      edges_buffer.push_back(b.x() + offset.x);
      edges_buffer.push_back(b.y() + offset.y);
      edges_buffer.push_back(b.z() + offset.z);

      //uvs [2]
       u = get(umap, opposite(halfedge(*he, *sm), *sm));
       v = get(vmap, opposite(halfedge(*he, *sm), *sm));

      edges_buffer.push_back(u);
      edges_buffer.push_back(v);
  }

  nb_face_verts = faces_buffer.size();
  nb_edge_verts = edges_buffer.size();
}

void Scene_textured_surface_mesh_item::common_constructor()
{
  cur_shading=FlatPlusEdges;
  is_selected=false;
  setTriangleContainer(0, new Tc(Vi::PROGRAM_WITH_TEXTURE,
                                 false));
  setEdgeContainer(1, new Ec(Three::mainViewer()->isOpenGL_4_3()
                             ? Vi::PROGRAM_SOLID_WIREFRAME
                             : Vi::PROGRAM_NO_SELECTION,
                             false));//bordures
  setEdgeContainer(0, new Ec(Vi::PROGRAM_WITH_TEXTURED_EDGES, false));//edges
  getTriangleContainer(0)->setTextureSize(QSize(d->texture.GetWidth(), d->texture.GetHeight()));
  getEdgeContainer(0)->setTextureSize(QSize(d->texture.GetWidth(), d->texture.GetHeight()));
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    initGL(viewer);
  }
  d->nb_face_verts = 0;
  d->nb_edge_verts = 0;
  d->nb_border_verts = 0;
  d->are_buffers_filled = false;
  d->are_border_filled = false;
}
Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item()
{
  d = new Scene_textured_surface_mesh_item_priv(this);
  common_constructor();

}

Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item(SMesh* const p)
{
  d = new Scene_textured_surface_mesh_item_priv(p,this);
  common_constructor();
}

Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item(const SMesh& p)
{
  d = new Scene_textured_surface_mesh_item_priv(p,this);
  common_constructor();
}

Scene_textured_surface_mesh_item::~Scene_textured_surface_mesh_item()
{
  delete d;
}

Scene_textured_surface_mesh_item*
Scene_textured_surface_mesh_item::clone() const {
  return new Scene_textured_surface_mesh_item(*d->sm);
}

// Load textured_polyhedron from .OFF file
bool
Scene_textured_surface_mesh_item::load(std::istream& in)
{
  std::cout<<"LOAD"<<std::endl;
  in >> *d->sm;
  invalidateOpenGLBuffers();
  return in && !isEmpty();
}

// Write textured_polyhedron to .OFF file
bool
Scene_textured_surface_mesh_item::save(std::ostream& out) const
{
  out << *d->sm;
  return (bool) out;
}

QString
Scene_textured_surface_mesh_item::toolTip() const
{
  if(!d->sm)
    return QString();

  return QObject::tr("<p>Textured polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
      .arg(this->name())
      .arg(num_vertices(*d->sm))
      .arg(num_halfedges(*d->sm)/2)
      .arg(num_faces(*d->sm))
      .arg(this->renderingModeName())
      .arg(this->color().name());
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_textured_surface_mesh_item::draw(CGAL::Three::Viewer_interface* viewer) const {

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
  getTriangleContainer(0)->draw(viewer, true);
}
void Scene_textured_surface_mesh_item::drawEdges(
    CGAL::Three::Viewer_interface* viewer) const {
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
  getEdgeContainer(0)->draw(viewer, true);

  if(viewer->isOpenGL_4_3())
  {

    QVector2D vp(viewer->width(), viewer->height());
    getEdgeContainer(1)->setViewport(vp);
    getEdgeContainer(1)->setWidth(4.0f);
  }
  getEdgeContainer(1)->setColor(QColor(Qt::blue));
  getEdgeContainer(1)->draw(viewer, true);
}


SMesh*
Scene_textured_surface_mesh_item::textured_face_graph()       { return d->sm; }
const SMesh*
Scene_textured_surface_mesh_item::textured_face_graph() const { return d->sm; }

bool
Scene_textured_surface_mesh_item::isEmpty() const {
  return (d->sm == 0) || d->sm->is_empty();
}

void
Scene_textured_surface_mesh_item::compute_bbox() const {
  const Point& p = d->sm->point(*vertices(*d->sm).begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());

  for(boost::graph_traits<SMesh>::vertex_iterator it = vertices(*d->sm).begin();
      it != vertices(*d->sm).end();
      ++it) {
    bbox = bbox + d->sm->point(*it).bbox();
  }
  setBbox(Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax()));
}
void
Scene_textured_surface_mesh_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  d->are_buffers_filled = false;
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  getEdgeContainer(0)->reset_vbos(ALL);
  getEdgeContainer(1)->reset_vbos(ALL);
}
void
Scene_textured_surface_mesh_item::selection_changed(bool p_is_selected)
{
  if(p_is_selected != is_selected)
  {
    is_selected = p_is_selected;
    Q_FOREACH(CGAL::QGLViewer*v, CGAL::QGLViewer::QGLViewerPool())
    {
      setBuffersInit(qobject_cast<Vi*>(v), false);
    }
    //to be replaced by a functor in the d-pointer when the merging is done
    if(p_is_selected)
      Q_EMIT selectionChanged();
  }
  else
    is_selected = p_is_selected;
}
void Scene_textured_surface_mesh_item::add_border_edges(std::vector<float> border_edges)
{
  d->border_edges_buffer = border_edges;
  d->are_border_filled = false;
  setBuffersFilled(false);
  itemChanged();
}

void Scene_textured_surface_mesh_item::initializeBuffers(Viewer_interface *v) const
{
    getTriangleContainer(0)->initializeBuffers(v);
    getEdgeContainer(0)->initializeBuffers(v);
    getTriangleContainer(0)->setFlatDataSize(d->nb_face_verts);
    getEdgeContainer(0)->setFlatDataSize(d->nb_edge_verts);
    d->edges_buffer.clear();
    d->edges_buffer.shrink_to_fit();
    d->faces_buffer.clear();
    d->faces_buffer.shrink_to_fit();

    getEdgeContainer(1)->initializeBuffers(v);
    getEdgeContainer(1)->setFlatDataSize(d->nb_border_verts);

    d->border_edges_buffer.clear();
    d->border_edges_buffer.shrink_to_fit();

}

void Scene_textured_surface_mesh_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!d->are_buffers_filled)
  {
    d->compute_normals_and_vertices();

    getTriangleContainer(0)->setTupleSize(8);
    getTriangleContainer(0)->getVbo(Tc::Flat_vertices)->stride = 8 * sizeof(float);
    getTriangleContainer(0)->allocate(
          Tc::Flat_vertices,
          d->faces_buffer.data(),
          static_cast<int>(d->faces_buffer.size()*sizeof(float)));

    getTriangleContainer(0)->getVbo(Tc::Flat_normals)->offset = 3 * sizeof(float);
    getTriangleContainer(0)->getVbo(Tc::Flat_normals)->stride = 8 * sizeof(float);
    getTriangleContainer(0)->allocate(
          Tc::Flat_normals,
          d->faces_buffer.data(),
          static_cast<int>(d->faces_buffer.size()*sizeof(float)));

    getTriangleContainer(0)->getVbo(Tc::Texture_map)->offset = 6 * sizeof(float);
    getTriangleContainer(0)->getVbo(Tc::Texture_map)->stride = 8 * sizeof(float);
    getTriangleContainer(0)->allocate(
          Tc::Texture_map,
          d->faces_buffer.data(),
          static_cast<int>(d->faces_buffer.size()*sizeof(float)));


    //Edges
    getEdgeContainer(0)->setTupleSize(5);
    getEdgeContainer(0)->getVbo(Ec::Vertices)->stride = 5 * sizeof(float);
    getEdgeContainer(0)->allocate(Ec::Vertices,
                                  d->edges_buffer.data(),
                                  static_cast<int>(d->edges_buffer.size()*sizeof(float)));
    getEdgeContainer(0)->getVbo(Ec::Texture_map)->offset =  3 * sizeof(float);
    getEdgeContainer(0)->getVbo(Ec::Texture_map)->stride = 5 * sizeof(float);
    getEdgeContainer(0)->allocate(Ec::Texture_map,
                                  d->edges_buffer.data(),
                                  static_cast<int>(d->edges_buffer.size()*sizeof(float)));

    getTriangleContainer(0)->getTexture()->setData(d->texture.GetData());
    getEdgeContainer(0)->getTexture()->setData(d->texture.GetData());
    d->are_buffers_filled = true;
  }
  if(!d->are_border_filled)
  {
    getEdgeContainer(1)->allocate(
          Ec::Vertices,
          d->border_edges_buffer.data(),
          static_cast<int>(d->border_edges_buffer.size()*sizeof(float)));
    d->nb_border_verts = d->border_edges_buffer.size();
    d->are_border_filled = true;
  }
  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}

