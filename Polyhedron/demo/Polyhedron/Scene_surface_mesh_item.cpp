#include "Scene_surface_mesh_item.h"

#include "Color_map.h"
#include <queue>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QApplication>
#include <QVariant>

//#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "triangulate_primitive.h"

typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::vertex_descriptor vertex_descriptor;
//Used to triangulate the AABB_Tree
class Primitive
{
public:
  // types
  typedef face_descriptor Id; // Id type
  typedef Scene_surface_mesh_item::Point Point; // point type
  typedef Scene_surface_mesh_item::Kernel::Triangle_3 Datum; // datum type

private:
  // member data
  Id m_it; // iterator
  Datum m_datum; // 3D triangle

  // constructor
public:
  Primitive() {}
  Primitive(Datum triangle, Id it)
    : m_it(it), m_datum(triangle)
  {
  }
public:
  Id& id() { return m_it; }
  const Id& id() const { return m_it; }
  Datum& datum() { return m_datum; }
  const Datum& datum() const { return m_datum; }

  /// Returns a point on the primitive
  Point reference_point() const { return m_datum.vertex(0); }
};


typedef CGAL::AABB_traits<Scene_surface_mesh_item::Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;

struct Scene_surface_mesh_item_priv{

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef std::vector<QColor> Color_vector;

  Scene_surface_mesh_item_priv(const Scene_surface_mesh_item& other, Scene_surface_mesh_item* parent):
    smesh_(new SMesh(*other.d->smesh_)),
    idx_data_(other.d->idx_data_),
    idx_edge_data_(other.d->idx_edge_data_)
  {
    item = parent;
  }

  Scene_surface_mesh_item_priv(SMesh* sm, Scene_surface_mesh_item *parent):
    smesh_(sm)
  {
    item = parent;
  }

  ~Scene_surface_mesh_item_priv()
  {
    delete smesh_;
  }

  void initialize_colors();
  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  void addFlatData(Point, Kernel::Vector_3, CGAL::Color *) const;
  void* get_aabb_tree();
  QList<Kernel::Triangle_3> triangulate_primitive(face_descriptor fit,
                                                  Kernel::Vector_3 normal);

  //! \brief triangulate_facet Triangulates a facet.
  //! \param fd a face_descriptor of the facet that needs to be triangulated.
  //! \param fnormals a property_map containing the normals of the mesh.
  //! \param fcolors a property_map containing the colors of the mesh
  //! \param im a property_map containing the indices of the vertices of the mesh
  //! \param index if true, the function will fill the index vector. If false, the function will
  //! fill the flat data vectors.
  void
  triangulate_facet(face_descriptor fd,
                    SMesh::Property_map<face_descriptor, Kernel::Vector_3 > *fnormals,
                    SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                    boost::property_map< SMesh, boost::vertex_index_t >::type* im,
                    bool index) const;
  void compute_elements();
  void checkFloat() const;

  enum VAOs {
   Flat_facets = 0,
   Smooth_facets,
   Edges,
   NbOfVaos
  };
  enum VBOs {
    Flat_vertices = 0,
    Smooth_vertices,
    Flat_normals,
    Smooth_normals,
    VColors,
    FColors,
    NbOfVbos
  };

  mutable bool has_fpatch_id;
  mutable bool floated;
  mutable bool has_vcolors;
  mutable bool has_fcolors;
  SMesh* smesh_;
  mutable bool is_filled;
  mutable bool isinit;
  mutable std::vector<unsigned int> idx_data_;
  std::vector<unsigned int> idx_edge_data_;
  mutable std::vector<cgal_gl_data> smooth_vertices;
  mutable std::vector<cgal_gl_data> smooth_normals;
  mutable std::vector<cgal_gl_data> flat_vertices;
  mutable std::vector<cgal_gl_data> flat_normals;
  mutable std::vector<cgal_gl_data> f_colors;
  mutable std::vector<cgal_gl_data> v_colors;
  mutable std::size_t nb_flat;
  mutable QOpenGLShaderProgram *program;
  Scene_surface_mesh_item *item;

  mutable SMesh::Property_map<face_descriptor,int> fpatch_id_map;
  mutable SMesh::Property_map<vertex_descriptor,int> v_selection_map;
  mutable SMesh::Property_map<face_descriptor,int> f_selection_map;

  Color_vector colors_;
  void computeElements() const;
};

const char* aabb_property_name = "Scene_surface_mesh_item aabb tree";

Scene_surface_mesh_item::Scene_surface_mesh_item()
  : CGAL::Three::Scene_item(Scene_surface_mesh_item_priv::NbOfVbos,Scene_surface_mesh_item_priv::NbOfVaos)
{
  d = new Scene_surface_mesh_item_priv(new SMesh(), this);
  d->floated = false;

  d->has_vcolors = false;
  d->has_fcolors = false;
  d->checkFloat();

  are_buffers_filled = false;
}

Scene_surface_mesh_item::Scene_surface_mesh_item(const Scene_surface_mesh_item& other)
  : CGAL::Three::Scene_item(Scene_surface_mesh_item_priv::NbOfVbos,Scene_surface_mesh_item_priv::NbOfVaos)
{
  d = new Scene_surface_mesh_item_priv(other, this);
  are_buffers_filled = false;
}

Scene_surface_mesh_item::Scene_surface_mesh_item(SMesh* sm)
  : CGAL::Three::Scene_item(Scene_surface_mesh_item_priv::NbOfVbos,Scene_surface_mesh_item_priv::NbOfVaos)
{
  d = new Scene_surface_mesh_item_priv(sm, this);
  d->floated = false;

  d->has_vcolors = false;
  d->has_fcolors = false;
  d->checkFloat();

  are_buffers_filled = false;
}

Scene_surface_mesh_item*
Scene_surface_mesh_item::clone() const
{ return new Scene_surface_mesh_item(*this); }


Scene_surface_mesh_item::Vertex_selection_map
Scene_surface_mesh_item::vertex_selection_map()
{
   if(! d->v_selection_map){
    d->v_selection_map = d->smesh_->add_property_map<vertex_descriptor,int>("v:selection").first;
   }
   return d->v_selection_map;
}

Scene_surface_mesh_item::Face_selection_map
Scene_surface_mesh_item::face_selection_map()
{
   if(! d->f_selection_map){
    d->f_selection_map = d->smesh_->add_property_map<face_descriptor,int>("f:selection").first;
   }
   return d->f_selection_map;
}

std::vector<QColor>&
Scene_surface_mesh_item::color_vector()
{
  return d->colors_;
}


void
Scene_surface_mesh_item::set_patch_id(SMesh::Face_index f,int i)const
{
  if(! d->has_fpatch_id){
    d->fpatch_id_map = d->smesh_->add_property_map<face_descriptor,int>("f:patch_id").first;
    d->has_fpatch_id = true;
  }
  d->fpatch_id_map[f] = i;
}

int
Scene_surface_mesh_item::patch_id(SMesh::Face_index f)const
{
  if(! d->has_fpatch_id){
    return -1;
  }
  return d->fpatch_id_map[f];
}

void Scene_surface_mesh_item_priv::addFlatData(Point p, Kernel::Vector_3 n, CGAL::Color *c) const
{
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();

  flat_vertices.push_back((cgal_gl_data)p.x()+offset[0]);
  flat_vertices.push_back((cgal_gl_data)p.y()+offset[1]);
  flat_vertices.push_back((cgal_gl_data)p.z()+offset[2]);

  flat_normals.push_back((cgal_gl_data)n.x());
  flat_normals.push_back((cgal_gl_data)n.y());
  flat_normals.push_back((cgal_gl_data)n.z());

  if(c != NULL)
  {
    f_colors.push_back((float)c->red()/255);
    f_colors.push_back((float)c->green()/255);
    f_colors.push_back((float)c->blue()/255);
  }
}

void Scene_surface_mesh_item_priv::compute_elements()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  smooth_vertices.clear();
  smooth_normals.clear();
  flat_vertices.clear();
  flat_normals.clear();
  f_colors.clear();
  v_colors.clear();
  idx_data_.clear();
  idx_data_.shrink_to_fit();

  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->add_property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, Kernel::Vector_3 > fnormals =
      smesh_->add_property_map<face_descriptor, Kernel::Vector_3 >("f:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*smesh_,fnormals);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  CGAL::Polygon_mesh_processing::compute_vertex_normals(*smesh_,vnormals);


  boost::property_map< SMesh, boost::vertex_index_t >::type
    im = get(boost::vertex_index, *smesh_);

  idx_data_.reserve(num_faces(*smesh_) * 3);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;


  BOOST_FOREACH(face_descriptor fd, faces(*smesh_))
  {
    if(is_triangle(halfedge(fd,*smesh_),*smesh_))
    {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
      {
        idx_data_.push_back(source(hd, *smesh_));
      }
    }
    else if(is_quad(halfedge(fd,*smesh_),*smesh_))
    {
      halfedge_descriptor hd = halfedge(fd,*smesh_);
      //1st half
        idx_data_.push_back(source(hd, *smesh_));
        idx_data_.push_back(source(next(hd, *smesh_), *smesh_));
        idx_data_.push_back(source(next(next(hd, *smesh_), *smesh_), *smesh_));

        //2nd half
        idx_data_.push_back(source(hd, *smesh_));
        idx_data_.push_back(source(next(next(hd, *smesh_), *smesh_), *smesh_));
        idx_data_.push_back(source(prev(hd, *smesh_), *smesh_));
    }
    else
    {
      triangulate_facet(fd, &fnormals, 0, &im, true);
    }
  }
  idx_edge_data_.clear();
  idx_edge_data_.shrink_to_fit();
  idx_edge_data_.reserve(num_edges(*smesh_) * 2);
  BOOST_FOREACH(edge_descriptor ed, edges(*smesh_))
  {
    idx_edge_data_.push_back(source(ed, *smesh_));
    idx_edge_data_.push_back(target(ed, *smesh_));
  }

  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();

  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
    smesh_->points();

  SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
    smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").first;

  SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
      smesh_->property_map<face_descriptor, CGAL::Color >("f:color").first;

  assert(positions.data() != NULL);
  assert(vnormals.data() != NULL);

  if(smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").second)
    has_vcolors = true;
  if(smesh_->property_map<face_descriptor, CGAL::Color >("f:color").second)
    has_fcolors = true;

  has_fpatch_id = smesh_->property_map<face_descriptor, int >("f:patch_id").second;

  if(has_fpatch_id && colors_.empty()){
    initialize_colors();
  }
//compute the Flat data
  flat_vertices.clear();
  flat_normals.clear();

  BOOST_FOREACH(face_descriptor fd, faces(*smesh_))
  {
    if(is_triangle(halfedge(fd,*smesh_),*smesh_))
    {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
      {
        Point p = positions[source(hd, *smesh_)];
        flat_vertices.push_back((cgal_gl_data)p.x()+offset.x);
        flat_vertices.push_back((cgal_gl_data)p.y()+offset.y);
        flat_vertices.push_back((cgal_gl_data)p.z()+offset.z);

        Kernel::Vector_3 n = fnormals[fd];
        flat_normals.push_back((cgal_gl_data)n.x());
        flat_normals.push_back((cgal_gl_data)n.y());
        flat_normals.push_back((cgal_gl_data)n.z());

        if(has_fpatch_id)
        {
          QColor c = item->color_vector()[fpatch_id_map[fd]];
          f_colors.push_back(c.redF());
          f_colors.push_back(c.greenF());
          f_colors.push_back(c.blueF());
        }
        else if(has_fcolors)
        {
          CGAL::Color c = fcolors[fd];
          f_colors.push_back((float)c.red()/255);
          f_colors.push_back((float)c.green()/255);
          f_colors.push_back((float)c.blue()/255);
        }
      }
    }
    else if(is_quad(halfedge(fd, *smesh_), *smesh_))
    {
      //1st half
      halfedge_descriptor hd = halfedge(fd, *smesh_);
      Point p = positions[source(hd, *smesh_)];
      Kernel::Vector_3 n = fnormals[fd];
      CGAL::Color *c;
      if(has_fcolors)
       c= &fcolors[fd];
      else
        c = 0;
      addFlatData(p,n,c);

      hd = next(halfedge(fd, *smesh_),*smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = next(next(halfedge(fd, *smesh_),*smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);
      //2nd half
      hd = halfedge(fd, *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = next(next(halfedge(fd, *smesh_),*smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = prev(halfedge(fd, *smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          , c);
    }
    else
    {
      triangulate_facet(fd, &fnormals, &fcolors, 0, false);
    }
  }

  if(has_vcolors)
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices(*smesh_))
    {
      CGAL::Color c = vcolors[vd];
      v_colors.push_back((float)c.red()/255);
      v_colors.push_back((float)c.green()/255);
      v_colors.push_back((float)c.blue()/255);
    }
  }

  if(floated)
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices(*smesh_))
    {
      Point p = positions[vd];
      smooth_vertices.push_back((cgal_gl_data)p.x()+offset.x);
      smooth_vertices.push_back((cgal_gl_data)p.y()+offset.y);
      smooth_vertices.push_back((cgal_gl_data)p.z()+offset.z);

      Kernel::Vector_3 n = vnormals[vd];
      smooth_normals.push_back((cgal_gl_data)n.x());
      smooth_normals.push_back((cgal_gl_data)n.y());
      smooth_normals.push_back((cgal_gl_data)n.z());

    }
  }
  QApplication::restoreOverrideCursor();
}


void Scene_surface_mesh_item_priv::initialize_colors()
{
  // Fill indices map and get max subdomain value
  int max = 0;
  int min = (std::numeric_limits<int>::max)();
  BOOST_FOREACH(face_descriptor fd, faces(*smesh_)){
    max = (std::max)(max, fpatch_id_map[fd]);
    min = (std::min)(min, fpatch_id_map[fd]);
  }

  colors_.clear();
  compute_color_map(item->color(), (std::max)(0, max + 1 - min),
                    std::back_inserter(colors_));
}

void Scene_surface_mesh_item_priv::initializeBuffers(CGAL::Three::Viewer_interface* viewer)const
{
  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
    smesh_->points();
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;
  //vao containing the data for the flat facets

  program = item->getShaderProgram(Scene_surface_mesh_item::PROGRAM_WITH_LIGHT, viewer);
  program->bind();

  item->vaos[Scene_surface_mesh_item_priv::Flat_facets]->bind();
  item->buffers[Scene_surface_mesh_item_priv::Flat_vertices].bind();
  item->buffers[Scene_surface_mesh_item_priv::Flat_vertices].allocate(flat_vertices.data(),
                             static_cast<int>(flat_vertices.size()*sizeof(cgal_gl_data)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",CGAL_GL_DATA,0,3);
  item->buffers[Scene_surface_mesh_item_priv::Flat_vertices].release();

  item->buffers[Scene_surface_mesh_item_priv::Flat_normals].bind();
  item->buffers[Scene_surface_mesh_item_priv::Flat_normals].allocate(flat_normals.data(),
                            static_cast<int>(flat_normals.size()*sizeof(cgal_gl_data)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals",CGAL_GL_DATA,0,3);
  item->buffers[Scene_surface_mesh_item_priv::Flat_normals].release();
  if(has_fcolors || has_fpatch_id)
  {
    item->buffers[Scene_surface_mesh_item_priv::FColors].bind();
    item->buffers[Scene_surface_mesh_item_priv::FColors].allocate(f_colors.data(),
                             static_cast<int>(f_colors.size()*sizeof(cgal_gl_data)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",CGAL_GL_DATA,0,3);
    item->buffers[Scene_surface_mesh_item_priv::FColors].release();
  }
  item->vaos[Scene_surface_mesh_item_priv::Flat_facets]->release();

  //vao containing the data for the smooth facets
  item->vaos[Scene_surface_mesh_item_priv::Smooth_facets]->bind();
  item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].bind();
  if(!(floated||viewer->offset() == qglviewer::Vec(0,0,0)))
  {
    item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].allocate(positions.data(),
                             static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  }
  else
  {
    item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].allocate(smooth_vertices.data(),
                               static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  }
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",CGAL_GL_DATA,0,3);
  item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].release();


  item->buffers[Scene_surface_mesh_item_priv::Smooth_normals].bind();
  if(!floated)
    item->buffers[Scene_surface_mesh_item_priv::Smooth_normals].allocate(vnormals.data(),
                                     static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  else
    item->buffers[Scene_surface_mesh_item_priv::Smooth_normals].allocate(smooth_normals.data(),
                              static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals",CGAL_GL_DATA,0,3);
  item->buffers[Scene_surface_mesh_item_priv::Smooth_normals].release();
  if(has_vcolors)
  {
    item->buffers[VColors].bind();
    item->buffers[VColors].allocate(v_colors.data(),
                             static_cast<int>(v_colors.size()*sizeof(cgal_gl_data)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",CGAL_GL_DATA,0,3);
    item->buffers[VColors].release();
  }
  item->vaos[Scene_surface_mesh_item_priv::Smooth_facets]->release();
  program->release();

  //vao for the edges
  program = item->getShaderProgram(Scene_surface_mesh_item::PROGRAM_WITHOUT_LIGHT, viewer);
  item->vaos[Scene_surface_mesh_item_priv::Edges]->bind();
  item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].bind();
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",CGAL_GL_DATA,0,3);
  item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].release();
  program->release();

  nb_flat = flat_vertices.size();
  smooth_vertices.resize(0);
  smooth_normals .resize(0);
  flat_vertices  .resize(0);
  flat_normals   .resize(0);
  f_colors       .resize(0);
  v_colors       .resize(0);
  smooth_vertices.shrink_to_fit();
  smooth_normals .shrink_to_fit();
  flat_vertices  .shrink_to_fit();
  flat_normals   .shrink_to_fit();
  f_colors       .shrink_to_fit();
  v_colors       .shrink_to_fit();

  item->are_buffers_filled = true;
}

void Scene_surface_mesh_item::draw(CGAL::Three::Viewer_interface *viewer) const
{
  glShadeModel(GL_SMOOTH);
  if(!are_buffers_filled)
  {
    d->compute_elements();
    d->initializeBuffers(viewer);
  }
  attribBuffers(viewer, PROGRAM_WITH_LIGHT);
  d->program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
  d->program->bind();

  if(renderingMode() == Gouraud)
  {
    vaos[Scene_surface_mesh_item_priv::Smooth_facets]->bind();
    if(is_selected)
      d->program->setAttributeValue("is_selected", true);
    else
      d->program->setAttributeValue("is_selected", false);
    if(!d->has_vcolors)
      d->program->setAttributeValue("colors", this->color());
    glDrawElements(GL_TRIANGLES, static_cast<GLuint>(d->idx_data_.size()),
                   GL_UNSIGNED_INT, d->idx_data_.data());
    vaos[Scene_surface_mesh_item_priv::Smooth_facets]->release();
  }
  else
  {
    vaos[Scene_surface_mesh_item_priv::Flat_facets]->bind();
    d->program->setAttributeValue("colors", this->color());
    if(is_selected)
      d->program->setAttributeValue("is_selected", true);
    else
      d->program->setAttributeValue("is_selected", false);
    if(!d->has_fcolors)
      d->program->setAttributeValue("colors", this->color());
    glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(d->nb_flat/3));
    vaos[Scene_surface_mesh_item_priv::Flat_facets]->release();
  }

  d->program->release();
}

void Scene_surface_mesh_item::drawEdges(CGAL::Three::Viewer_interface *viewer) const
{
  if(!are_buffers_filled)
  {
    d->compute_elements();
    d->initializeBuffers(viewer);
  }
 attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
 d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
 d->program->bind();
 vaos[Scene_surface_mesh_item_priv::Edges]->bind();
 d->program->setAttributeValue("colors", QColor(0,0,0));
 if(is_selected)
   d->program->setAttributeValue("is_selected", true);
 else
   d->program->setAttributeValue("is_selected", false);
 glDrawElements(GL_LINES, static_cast<GLuint>(d->idx_edge_data_.size()),
                GL_UNSIGNED_INT, d->idx_edge_data_.data());
 vaos[Scene_surface_mesh_item_priv::Edges]->release();
 d->program->release();
}

void Scene_surface_mesh_item::drawPoints(CGAL::Three::Viewer_interface *viewer) const
{
  if(!are_buffers_filled)
  {
    d->compute_elements();
    d->initializeBuffers(viewer);
  }
 attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
 d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
 d->program->bind();
 vaos[Scene_surface_mesh_item_priv::Edges]->bind();
 d->program->setAttributeValue("colors", QColor(0,0,0));
 if(is_selected)
   d->program->setAttributeValue("is_selected", true);
 else
   d->program->setAttributeValue("is_selected", false);
 glDrawElements(GL_POINTS, static_cast<GLuint>(d->idx_edge_data_.size()),
                GL_UNSIGNED_INT, d->idx_edge_data_.data());
 vaos[Scene_surface_mesh_item_priv::Edges]->release();
 d->program->release();
}

void
Scene_surface_mesh_item::selection_changed(bool p_is_selected)
{
  if(p_is_selected != is_selected)
  {
    is_selected = p_is_selected;
  }
}

bool
Scene_surface_mesh_item::supportsRenderingMode(RenderingMode m) const
{ return (m == FlatPlusEdges || m == Wireframe || m == Flat || m == Gouraud || m == Points); }

CGAL::Three::Scene_item::Bbox Scene_surface_mesh_item::bbox() const
{
 if(!is_bbox_computed)
   compute_bbox();
 return _bbox;
}

bool
Scene_surface_mesh_item::isEmpty() const
{

  return num_vertices(*d->smesh_)==0;
}

QString Scene_surface_mesh_item::toolTip() const
{
  return QObject::tr("<p>Surface_mesh <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of faces: %4</p>")
    .arg(this->name())
    .arg(num_vertices(*d->smesh_))
    .arg(num_edges(*d->smesh_))
    .arg(num_faces(*d->smesh_))
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

void Scene_surface_mesh_item_priv::checkFloat()const
{
#if CGAL_IS_FLOAT == 1
  floated = true;
#endif
}

void
Scene_surface_mesh_item_priv::triangulate_facet(face_descriptor fd,
                                           SMesh::Property_map<face_descriptor, Kernel::Vector_3> *fnormals,
                                           SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                                           boost::property_map< SMesh, boost::vertex_index_t >::type *im,
                                           bool index) const
{
  //Computes the normal of the facet
  Kernel::Vector_3 normal = get(*fnormals, fd);

  //check if normal contains NaN values
  if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
  {
    qDebug()<<"Warning : normal is not valid. Facet not displayed";
    return;
  }

  typedef FacetTriangulator<SMesh, Kernel, boost::graph_traits<SMesh>::vertex_descriptor> FT;
  double diagonal;
  if(item->diagonalBbox() != std::numeric_limits<double>::infinity())
    diagonal = item->diagonalBbox();
  else
    diagonal = 0.0;
  FT triangulation(fd,normal,smesh_,diagonal);
  //iterates on the internal faces
  for(FT::CDT::Finite_faces_iterator
       ffit = triangulation.cdt->finite_faces_begin(),
       end = triangulation.cdt->finite_faces_end();
       ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;
    //add the vertices to the positions
      //adds the vertices, normals and colors to the appropriate vectors
    if(!index)
    {
      CGAL::Color* color;
      if(has_fcolors)
        color = &(*fcolors)[fd];
      else
        color = 0;

      addFlatData(ffit->vertex(0)->point(),
                  (*fnormals)[fd],
                  color);
      addFlatData(ffit->vertex(1)->point(),
                  (*fnormals)[fd],
                  color);

      addFlatData(ffit->vertex(2)->point(),
                  (*fnormals)[fd],
                  color);
    }
    //adds the indices to the appropriate vector
    else
    {
      idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(0)]]);
      idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(1)]]);
      idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(2)]]);
    }

  }
}
void delete_aabb_tree(Scene_surface_mesh_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
        if(tree) {
            delete tree;
            tree = 0;
        }
        item->setProperty(aabb_property_name, QVariant());
    }
}

Scene_surface_mesh_item::~Scene_surface_mesh_item()
{
  delete_aabb_tree(this);
  delete d;
}
Scene_surface_mesh_item::SMesh* Scene_surface_mesh_item::polyhedron() { return d->smesh_; }
const Scene_surface_mesh_item::SMesh* Scene_surface_mesh_item::polyhedron() const { return d->smesh_; }

void Scene_surface_mesh_item::compute_bbox()const
{
  SMesh::Property_map<vertex_descriptor, Point> pprop = d->smesh_->points();
  CGAL::Bbox_3 bbox;

  BOOST_FOREACH(vertex_descriptor vd,vertices(*d->smesh_))
  {
    bbox = bbox + pprop[vd].bbox();
  }
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax());

}

void* Scene_surface_mesh_item_priv::get_aabb_tree()
{
  QVariant aabb_tree_property = item->property(aabb_property_name);
  if(aabb_tree_property.isValid()) {
    void* ptr = aabb_tree_property.value<void*>();
    return static_cast<Input_facets_AABB_tree*>(ptr);
  }
  else {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    SMesh* sm = item->polyhedron();
    if(sm) {

      Input_facets_AABB_tree* tree =
          new Input_facets_AABB_tree();
      int index =0;
      BOOST_FOREACH( face_descriptor f, faces(*sm))
      {
        //if face not triangle, triangulate corresponding primitive before adding it to the tree
        if(!CGAL::is_triangle(halfedge(f, *sm), *sm))
        {
          Kernel::Plane_3 face_plane = Kernel::Plane_3(sm->point(target(halfedge(f, *sm), *sm)),
                                                       sm->point(target(next(halfedge(f, *sm), *sm), *sm)),
                                                       sm->point(target(next(next(halfedge(f, *sm), *sm), *sm), *sm)));
          Kernel::Vector_3 normal = face_plane.orthogonal_vector(); //initialized in compute_normals_and_vertices
          index +=3;
          Q_FOREACH(Kernel::Triangle_3 triangle, triangulate_primitive(f,normal))
          {
            Primitive primitive(triangle, f);
            tree->insert(primitive);
          }
        }
        else
        {
          Kernel::Triangle_3 triangle(
                sm->point(target(halfedge(f, *sm), *sm)),
                sm->point(target(next(halfedge(f, *sm), *sm), *sm)),
                sm->point(target(next(next(halfedge(f, *sm), *sm), *sm), *sm))
                );
          Primitive primitive(triangle, f);
          tree->insert(primitive);
        }
      }
      item->setProperty(aabb_property_name,
                        QVariant::fromValue<void*>(tree));
      QApplication::restoreOverrideCursor();
      return tree;
    }
    else return 0;
  }
}


void
Scene_surface_mesh_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
  SMesh *sm = d->smesh_;
  std::size_t vertex_to_emit = 0;
  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Intersection_and_primitive_id<Kernel::Ray_3>::Type Object_and_primitive_id;

  Tree* aabb_tree = static_cast<Tree*>(d->get_aabb_tree());
  if(aabb_tree)
  {
    const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
    const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
    const Kernel::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Object_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));
    Intersections::iterator closest = intersections.begin();
    if(closest != intersections.end())
    {

      const Kernel::Point_3* closest_point =
          boost::get<Kernel::Point_3>(&(closest->first));
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const Kernel::Point_3* it_point =
              boost::get<Kernel::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point) {
        face_descriptor selected_face = closest->second;

        // The computation of the nearest vertex may be costly.  Only
        // do it if some objects are connected to the signal
        // 'selected_vertex'.
        if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
        {

          SMesh::Halfedge_around_face_circulator he_it(sm->halfedge(selected_face),*sm), around_end(he_it);

          vertex_descriptor v = sm->target(*he_it), nearest_v = v;

          Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                      sm->point(v));
          while(++he_it != around_end) {
            v = sm->target(*he_it);
            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                            sm->point(v));
            if(new_sq_dist < sq_dist) {
              sq_dist = new_sq_dist;
              nearest_v = v;
            }
          }
          //bottleneck
          vertex_to_emit = static_cast<std::size_t>(nearest_v);
        }

        if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
           || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
        {
          SMesh::Halfedge_around_face_circulator he_it(sm->halfedge(selected_face),*sm), around_end(he_it);

          halfedge_descriptor nearest_h = *he_it;
          Kernel::FT sq_dist =
              CGAL::squared_distance(*closest_point,
                                     Kernel::Segment_3(sm->point(sm->target(*he_it)),
                                                       sm->point(
                                                         sm->target(
                                                           sm->opposite(*he_it)))));

          while(++he_it != around_end)
          {
            Kernel::FT new_sq_dist =
                CGAL::squared_distance(*closest_point,
                                       Kernel::Segment_3(sm->point(sm->target(*he_it)),
                                                         sm->point(
                                                           sm->target(
                                                             sm->opposite(*he_it)))));
            if(new_sq_dist < sq_dist) {
              sq_dist = new_sq_dist;
              nearest_h = *he_it;
            }
          }
          std::size_t s_nearest_h = static_cast<std::size_t>(nearest_h);
          std::size_t s_nearest_e = static_cast<std::size_t>(nearest_h)/2;
          Q_EMIT selected_halfedge(reinterpret_cast<void*>(s_nearest_h));
          Q_EMIT selected_edge(reinterpret_cast<void*>(s_nearest_e));
        }
        Q_EMIT selected_vertex(reinterpret_cast<void*>(vertex_to_emit));
        std::size_t s_selected_f = static_cast<std::size_t>(selected_face);
        Q_EMIT selected_facet(reinterpret_cast<void*>(s_selected_f));
      }
    }
  }
  Scene_item::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
  Q_EMIT selection_done();
}

void Scene_surface_mesh_item::invalidateOpenGLBuffers()
{
  Q_EMIT item_is_about_to_be_changed();
  delete_aabb_tree(this);
  d->smesh_->collect_garbage();
  are_buffers_filled = false;
}


QList<Scene_surface_mesh_item::Kernel::Triangle_3> Scene_surface_mesh_item_priv::triangulate_primitive(face_descriptor fit,
                                                Scene_surface_mesh_item::Kernel::Vector_3 normal)
{
  typedef FacetTriangulator<SMesh, Kernel, boost::graph_traits<SMesh>::vertex_descriptor> FT;
  //The output list
  QList<Kernel::Triangle_3> res;
  //check if normal contains NaN values
  if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
  {
    qDebug()<<"Warning in triangulation of the selection item: normal contains NaN values and is not valid.";
    return QList<Kernel::Triangle_3>();
  }
  double diagonal;
  if(item->diagonalBbox() != std::numeric_limits<double>::infinity())
    diagonal = item->diagonalBbox();
  else
    diagonal = 0.0;
  FT triangulation(fit,normal,smesh_,diagonal);
  //iterates on the internal faces to add the vertices to the positions
  //and the normals to the appropriate vectors
  for( FT::CDT::Finite_faces_iterator
      ffit = triangulation.cdt->finite_faces_begin(),
      end = triangulation.cdt->finite_faces_end();
      ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;


    res << Kernel::Triangle_3(ffit->vertex(0)->point(),
                              ffit->vertex(1)->point(),
                              ffit->vertex(2)->point());

  }
  return res;
}

void Scene_surface_mesh_item::invalidate_aabb_tree()
{
 delete_aabb_tree(this);
}


bool Scene_surface_mesh_item::intersect_face(double orig_x,
                                           double orig_y,
                                           double orig_z,
                                           double dir_x,
                                           double dir_y,
                                           double dir_z,
                                           const face_descriptor &f)
{
  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Object_and_primitive_id Object_and_primitive_id;

  Tree* aabb_tree = static_cast<Tree*>(d->get_aabb_tree());
  if(aabb_tree)
  {
    const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
    const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
    const Kernel::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Object_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));
    Intersections::iterator closest = intersections.begin();
    if(closest != intersections.end())
    {
      const Kernel::Point_3* closest_point =
          CGAL::object_cast<Kernel::Point_3>(&closest->first);
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const Kernel::Point_3* it_point =
              CGAL::object_cast<Kernel::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point)
      {
        face_descriptor intersected_face = closest->second;
        return intersected_face == f;
      }
    }
  }
  return false;

}
