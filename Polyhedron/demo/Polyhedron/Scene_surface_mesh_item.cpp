#include "Scene_surface_mesh_item.h"

#include <queue>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QApplication>

#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "triangulate_primitive.h"

typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::vertex_descriptor vertex_descriptor;


struct Scene_surface_mesh_item_priv{

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;


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

  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  void addFlatData(Point, Kernel::Vector_3, CGAL::Color *) const;

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
  mutable QOpenGLShaderProgram *program;
  Scene_surface_mesh_item *item;

};
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
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    d->smesh_->add_property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, Kernel::Vector_3 > fnormals =
      d->smesh_->add_property_map<face_descriptor, Kernel::Vector_3 >("v:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*d->smesh_,fnormals);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  CGAL::Polygon_mesh_processing::compute_vertex_normals(*d->smesh_,vnormals);


  boost::property_map< SMesh, boost::vertex_index_t >::type
    im = get(boost::vertex_index, *d->smesh_);

  d->idx_data_.reserve(num_faces(*d->smesh_) * 3);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;



  BOOST_FOREACH(face_descriptor fd, faces(*d->smesh_))
  {
    if(is_triangle(halfedge(fd,*d->smesh_),*d->smesh_))
    {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *d->smesh_),*d->smesh_))
      {
        d->idx_data_.push_back(im[source(hd, *d->smesh_)]);
      }
    }
    else if(is_quad(halfedge(fd,*d->smesh_),*d->smesh_))
    {
      halfedge_descriptor hd = halfedge(fd,*d->smesh_);
      //1st half
        d->idx_data_.push_back(im[source(hd, *d->smesh_)]);
        d->idx_data_.push_back(im[source(next(hd, *d->smesh_), *d->smesh_)]);
        d->idx_data_.push_back(im[source(next(next(hd, *d->smesh_), *d->smesh_), *d->smesh_)]);

        //2nd half
        d->idx_data_.push_back(im[source(hd, *d->smesh_)]);
        d->idx_data_.push_back(im[source(next(next(hd, *d->smesh_), *d->smesh_), *d->smesh_)]);
        d->idx_data_.push_back(im[source(prev(hd, *d->smesh_), *d->smesh_)]);
    }
    else
    {
      d->triangulate_facet(fd, &fnormals, 0, &im, true);
    }
  }

  d->idx_edge_data_.reserve(num_edges(*d->smesh_) * 2);
  BOOST_FOREACH(edge_descriptor ed, edges(*d->smesh_))
  {
    d->idx_edge_data_.push_back(im[source(ed, *d->smesh_)]);
    d->idx_edge_data_.push_back(im[target(ed, *d->smesh_)]);
  }

  d->compute_elements();
  are_buffers_filled = false;
}

Scene_surface_mesh_item*
Scene_surface_mesh_item::clone() const
{ return new Scene_surface_mesh_item(*this); }

void Scene_surface_mesh_item_priv::addFlatData(Point p, Kernel::Vector_3 n, CGAL::Color *c) const
{

  flat_vertices.push_back((cgal_gl_data)p.x());
  flat_vertices.push_back((cgal_gl_data)p.y());
  flat_vertices.push_back((cgal_gl_data)p.z());

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
  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
    smesh_->points();
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, Kernel::Vector_3 > fnormals =
      smesh_->add_property_map<face_descriptor, Kernel::Vector_3 >("v:normal").first;

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
        flat_vertices.push_back((cgal_gl_data)p.x());
        flat_vertices.push_back((cgal_gl_data)p.y());
        flat_vertices.push_back((cgal_gl_data)p.z());

        Kernel::Vector_3 n = fnormals[fd];
        flat_normals.push_back((cgal_gl_data)n.x());
        flat_normals.push_back((cgal_gl_data)n.y());
        flat_normals.push_back((cgal_gl_data)n.z());

        if(has_fcolors)
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
      smooth_vertices.push_back((cgal_gl_data)p.x());
      smooth_vertices.push_back((cgal_gl_data)p.y());
      smooth_vertices.push_back((cgal_gl_data)p.z());

      Kernel::Vector_3 n = vnormals[vd];
      smooth_normals.push_back((cgal_gl_data)n.x());
      smooth_normals.push_back((cgal_gl_data)n.y());
      smooth_normals.push_back((cgal_gl_data)n.z());

    }
  }
  QApplication::restoreOverrideCursor();
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
  if(has_fcolors)
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
  if(!floated)
    item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].allocate(positions.data(),
                             static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  else
    item->buffers[Scene_surface_mesh_item_priv::Smooth_vertices].allocate(smooth_vertices.data(),
                               static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
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
  item->are_buffers_filled = true;
}

void Scene_surface_mesh_item::draw(CGAL::Three::Viewer_interface *viewer) const
{
  glShadeModel(GL_SMOOTH);
  if(!are_buffers_filled)
    d->initializeBuffers(viewer);
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
    glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(d->flat_vertices.size()/3));
    vaos[Scene_surface_mesh_item_priv::Flat_facets]->release();
  }

  d->program->release();
}

void Scene_surface_mesh_item::drawEdges(CGAL::Three::Viewer_interface *viewer) const
{
 if(!are_buffers_filled)
   d->initializeBuffers(viewer);
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

void Scene_surface_mesh_item::drawPoints(CGAL::Three::Viewer_interface *) const
{

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
{ return (m == FlatPlusEdges || m == Wireframe || m == Flat || m == Gouraud); }

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

Scene_surface_mesh_item::~Scene_surface_mesh_item()
{
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

