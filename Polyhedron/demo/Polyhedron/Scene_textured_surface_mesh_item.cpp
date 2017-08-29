#include "Scene_textured_surface_mesh_item.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <QApplication>
#include <QObject>

typedef EPICK  ::Point_3 Point;

struct Scene_textured_surface_mesh_item_priv
{
  Scene_textured_surface_mesh_item_priv(Scene_textured_surface_mesh_item* parent)
    :sm(new SMesh), textureId(-1)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    uv = sm->add_property_map<halfedge_descriptor,std::pair<float, float> >("h:uv",std::make_pair(0.0f,0.0f)).first;

  }
  Scene_textured_surface_mesh_item_priv(const SMesh& p, Scene_textured_surface_mesh_item* parent)
    : sm(new SMesh(p)),textureId(-1),smooth_shading(true)

  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    uv = sm->add_property_map<halfedge_descriptor,std::pair<float, float> >("h:uv",std::make_pair(0.0f,0.0f)).first;
  }
  Scene_textured_surface_mesh_item_priv(SMesh* const p,Scene_textured_surface_mesh_item* parent)
    :sm(p),textureId(-1),smooth_shading(true)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    uv = sm->add_property_map<halfedge_descriptor,std::pair<float, float> >("h:uv",std::make_pair(0.0f,0.0f)).first;
  }

  ~Scene_textured_surface_mesh_item_priv()
  {
    delete sm;
  }

  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices(void) const;

  enum VAOs {
    Facets=0,
    Edges,
    Border_edges,
    NbOfVaos
  };
  enum VBOs {
    B_Facets=0,
    B_Edges,
    B_Border_Edges,
    NbOfVbos
  };

  SMesh* sm;
  Texture texture;
  SMesh::Property_map<halfedge_descriptor,std::pair<float, float> > uv;
  mutable GLuint textureId;
  mutable QOpenGLShaderProgram* program;
  //[Px][Py][Pz][Nx][Ny][Nz][u][v]
  mutable std::vector<float> faces_buffer;
  //[Px][Py][Pz][u][v]
  mutable std::vector<float> edges_buffer;
  //[Px][Py][Pz]
  mutable std::vector<float> border_edges_buffer;
  bool smooth_shading;

  Scene_textured_surface_mesh_item* item;
  typedef Scene_textured_surface_mesh_item I;
};
void Scene_textured_surface_mesh_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer = 0) const
{
  if(GLuint(-1) == textureId) {
    viewer->glGenTextures(1, &textureId);
  }

  //Faces
  program = item->getShaderProgram(I::PROGRAM_WITH_TEXTURE, viewer);
  program->bind();
  item->vaos[Facets]->bind();
  item->buffers[B_Facets].bind();
  item->buffers[B_Facets].allocate(faces_buffer.data(),
                                   static_cast<int>(faces_buffer.size()*sizeof(float)));

  program->enableAttributeArray("vertex");
  program->enableAttributeArray("normal");
  program->enableAttributeArray("v_texCoord");

  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3, 8 * sizeof(float));
  program->setAttributeBuffer("normal", GL_FLOAT, 3 * sizeof(float), 3, 8 * sizeof(float));
  program->setAttributeBuffer("v_texCoord", GL_FLOAT, 6 * sizeof(float), 2, 8 * sizeof(float));

  item->buffers[B_Facets].release();
  item->vaos[Facets]->release();
  program->release();

  //Edges
  program = item->getShaderProgram(I::PROGRAM_WITH_TEXTURED_EDGES, viewer);
  program->bind();
  item->vaos[Edges]->bind();
  item->buffers[B_Edges].bind();
  item->buffers[B_Edges].allocate(edges_buffer.data(),
                                   static_cast<int>(edges_buffer.size()*sizeof(float)));

  program->enableAttributeArray("vertex");
  program->enableAttributeArray("v_texCoord");

  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3, 5 * sizeof(float));
  program->setAttributeBuffer("v_texCoord", GL_FLOAT, 3 * sizeof(float), 2, 5 * sizeof(float));

  item->buffers[B_Edges].release();
  item->vaos[Edges]->release();
  program->release();




  viewer->glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  viewer->glActiveTexture(GL_TEXTURE0);
  viewer->glBindTexture(GL_TEXTURE_2D, textureId);
  viewer->glTexImage2D(GL_TEXTURE_2D,
                       0,
                       GL_RGB,
                       texture.GetWidth(),
                       texture.GetHeight(),
                       0,
                       GL_RGB,
                       GL_UNSIGNED_BYTE,
                       texture.GetData());
  viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  item->are_buffers_filled = true;
}

void
Scene_textured_surface_mesh_item_priv::compute_normals_and_vertices(void) const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  faces_buffer.resize(0);


  typedef boost::graph_traits<SMesh>::face_iterator face_iterator;
  typedef boost::graph_traits<SMesh>::face_iterator face_iterator;
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();

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
      const float u = get(uv, *he).first;
      const float v = get(uv, *he).second;
      faces_buffer.push_back(u);
      faces_buffer.push_back(v);
    }
  }

  //Edges
  typedef EPICK::Point_3		Point;
  typedef SMesh::Edge_iterator	Edge_iterator;

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
      float u = get(uv, halfedge(*he, *sm)).first;
      float v = get(uv, halfedge(*he, *sm)).second;

      edges_buffer.push_back(u);
      edges_buffer.push_back(v);
      //position [3]
      edges_buffer.push_back(b.x() + offset.x);
      edges_buffer.push_back(b.y() + offset.y);
      edges_buffer.push_back(b.z() + offset.z);

      //uvs [2]
       u = get(uv, opposite(halfedge(*he, *sm), *sm)).first;
       v = get(uv, opposite(halfedge(*he, *sm), *sm)).second;

      edges_buffer.push_back(u);
      edges_buffer.push_back(v);

  }


  QApplication::restoreOverrideCursor();
}

Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item()
  : Scene_item(Scene_textured_surface_mesh_item_priv::NbOfVbos,Scene_textured_surface_mesh_item_priv::NbOfVaos)
{
  cur_shading=FlatPlusEdges;
  is_selected=false;
  d = new Scene_textured_surface_mesh_item_priv(this);
  invalidateOpenGLBuffers();
}

Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item(SMesh* const p)
  : Scene_item(Scene_textured_surface_mesh_item_priv::NbOfVbos,Scene_textured_surface_mesh_item_priv::NbOfVaos)
{
  cur_shading=FlatPlusEdges;
  is_selected=false;
  d = new Scene_textured_surface_mesh_item_priv(p,this);
  invalidateOpenGLBuffers();
}

Scene_textured_surface_mesh_item::Scene_textured_surface_mesh_item(const SMesh& p)
  : Scene_item(Scene_textured_surface_mesh_item_priv::NbOfVbos,Scene_textured_surface_mesh_item_priv::NbOfVaos)
{
  cur_shading=FlatPlusEdges;
  is_selected=false;
  d = new Scene_textured_surface_mesh_item_priv(p,this);
  invalidateOpenGLBuffers();
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

  if(!are_buffers_filled)
  {
    d->compute_normals_and_vertices();
    d->initializeBuffers(viewer);
  }
  vaos[Scene_textured_surface_mesh_item_priv::Facets]->bind();
  viewer->glActiveTexture(GL_TEXTURE0);
  viewer->glBindTexture(GL_TEXTURE_2D, d->textureId);
  attribBuffers(viewer, PROGRAM_WITH_TEXTURE);
  d->program=getShaderProgram(PROGRAM_WITH_TEXTURE);
  d->program->bind();
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->faces_buffer.size()/8));
  //Clean-up
  d->program->release();
  vaos[Scene_textured_surface_mesh_item_priv::Facets]->release();

}
void Scene_textured_surface_mesh_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
  if(!are_buffers_filled)
    d->initializeBuffers(viewer);

  vaos[Scene_textured_surface_mesh_item_priv::Edges]->bind();
  viewer->glActiveTexture(GL_TEXTURE0);
  viewer->glBindTexture(GL_TEXTURE_2D, d->textureId);
  attribBuffers(viewer, PROGRAM_WITH_TEXTURED_EDGES);

  d->program=getShaderProgram(PROGRAM_WITH_TEXTURED_EDGES);
  d->program->bind();
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->edges_buffer.size()/5));

  vaos[Scene_textured_surface_mesh_item_priv::Edges]->release();
  d->program->release();

  vaos[Scene_textured_surface_mesh_item_priv::Border_edges]->bind();
  attribBuffers(viewer, PROGRAM_NO_SELECTION);
  d->program=getShaderProgram(PROGRAM_NO_SELECTION);
  d->program->bind();
  viewer->glLineWidth(4.0);
  d->program->setAttributeValue("colors", QColor(Qt::blue));
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->border_edges_buffer.size()/3));
  viewer->glLineWidth(1.0);
  //Clean-up
  d->program->release();
  vaos[Scene_textured_surface_mesh_item_priv::Border_edges]->release();
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
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax());
}
void
Scene_textured_surface_mesh_item::invalidateOpenGLBuffers()
{
  are_buffers_filled = false;
  compute_bbox();
}
void
Scene_textured_surface_mesh_item::selection_changed(bool p_is_selected)
{
  if(p_is_selected != is_selected)
  {
    is_selected = p_is_selected;
    initializeBuffers();
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
  d->program=getShaderProgram(PROGRAM_NO_SELECTION);
  d->program->bind();
  vaos[Scene_textured_surface_mesh_item_priv::Border_edges]->bind();
  buffers[Scene_textured_surface_mesh_item_priv::B_Border_Edges].bind();
  buffers[Scene_textured_surface_mesh_item_priv::B_Border_Edges].allocate(d->border_edges_buffer.data(),
                      static_cast<int>(d->border_edges_buffer.size()*sizeof(float)));
  d->program->enableAttributeArray("vertex");
  d->program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  d->program->disableAttributeArray("colors");
  buffers[Scene_textured_surface_mesh_item_priv::B_Border_Edges].release();
  vaos[Scene_textured_surface_mesh_item_priv::Border_edges]->release();

  d->program->release();
  itemChanged();

}
