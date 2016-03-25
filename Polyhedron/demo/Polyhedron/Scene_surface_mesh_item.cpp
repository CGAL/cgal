#include "Scene_surface_mesh_item.h"

#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

Scene_surface_mesh_item::Scene_surface_mesh_item(const Scene_surface_mesh_item& other)
  : CGAL::Three::Scene_item(NbOfVbos,NbOfVaos),
    smesh_(new SMesh(*other.smesh_)),
    idx_data_(other.idx_data_),
    idx_edge_data_(other.idx_edge_data_)
{
}

Scene_surface_mesh_item::Scene_surface_mesh_item(SMesh* sm)
  : CGAL::Three::Scene_item(NbOfVbos,NbOfVaos),
    smesh_(sm)
{
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->add_property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  CGAL::Polygon_mesh_processing::compute_vertex_normals(*smesh_,vnormals);


  boost::property_map< SMesh, boost::vertex_index_t >::type
    im = boost::get(boost::vertex_index, *smesh_);

  idx_data_.reserve(num_faces(*smesh_) * 3);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;



  BOOST_FOREACH(face_descriptor fd, faces(*smesh_))
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
    {
      idx_data_.push_back(im[source(hd, *smesh_)]);
    }
  }

  idx_edge_data_.reserve(num_edges(*smesh_) * 2);
  BOOST_FOREACH(edge_descriptor ed, edges(*smesh_))
  {
    idx_edge_data_.push_back(im[source(ed, *smesh_)]);
    idx_edge_data_.push_back(im[target(ed, *smesh_)]);
  }
}

Scene_surface_mesh_item*
Scene_surface_mesh_item::clone() const
{ return new Scene_surface_mesh_item(*this); }


void Scene_surface_mesh_item::initializeBuffers(CGAL::Three::Viewer_interface* viewer)const
{
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
    smesh_->points();
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, Kernel::Vector_3 > fnormals =
      smesh_->add_property_map<face_descriptor, Kernel::Vector_3 >("v:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*smesh_,fnormals);
  assert(positions.data() != NULL);
  assert(vnormals.data() != NULL);

//compute the Flat data
  flat_vertices.clear();
  flat_normals.clear();
  BOOST_FOREACH(face_descriptor fd, faces(*smesh_))
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
    {
      Point p = positions[source(hd, *smesh_)];
      flat_vertices.push_back((gl_data)p.x());
      flat_vertices.push_back((gl_data)p.y());
      flat_vertices.push_back((gl_data)p.z());
      Kernel::Vector_3 n = fnormals[fd];
      flat_normals.push_back((gl_data)n.x());
      flat_normals.push_back((gl_data)n.y());
      flat_normals.push_back((gl_data)n.z());
    }
  }


  //vao containing the data for the flat facets

  program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
  program->bind();

  vaos[Flat_facets]->bind();
  buffers[Flat_vertices].bind();
  buffers[Flat_vertices].allocate(flat_vertices.data(),
                             static_cast<int>(flat_vertices.size()*sizeof(gl_data)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_DATA,0,3);
  buffers[Flat_vertices].release();



  buffers[Flat_normals].bind();
  buffers[Flat_normals].allocate(flat_normals.data(),
                            static_cast<int>(flat_normals.size()*sizeof(gl_data)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals",GL_DATA,0,3);
  buffers[Flat_normals].release();
  vaos[Flat_facets]->release();

  //vao containing the data for the smooth facets

  vaos[Smooth_facets]->bind();
  buffers[Smooth_vertices].bind();
  buffers[Smooth_vertices].allocate(positions.data(),
                             static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_DATA,0,3);
  buffers[Smooth_vertices].release();



  buffers[Smooth_normals].bind();
  buffers[Smooth_normals].allocate(vnormals.data(),
                            static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals",GL_DATA,0,3);
  buffers[Smooth_normals].release();
  vaos[Smooth_facets]->release();
  program->release();

  //vao for the edges
  program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
  vaos[Edges]->bind();
  buffers[Smooth_vertices].bind();
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_DATA,0,3);
  buffers[Smooth_vertices].release();
  program->release();
  are_buffers_filled = true;
}

void Scene_surface_mesh_item::draw(CGAL::Three::Viewer_interface *viewer) const
{

  if(!are_buffers_filled)
    initializeBuffers(viewer);
  attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
  program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
  program->bind();

  if(renderingMode() == Gouraud)
  {
    vaos[Smooth_facets]->bind();
    program->setAttributeValue("colors", this->color());
    glDrawElements(GL_TRIANGLES, idx_data_.size(),
                   GL_UNSIGNED_INT, idx_data_.data());
    vaos[Smooth_facets]->release();
  }
  else
  {
    vaos[Flat_facets]->bind();
    program->setAttributeValue("colors", this->color());
    glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(flat_vertices.size()/3));
    vaos[Flat_facets]->release();
  }

  program->release();
}

void Scene_surface_mesh_item::draw_edges(CGAL::Three::Viewer_interface *viewer) const
{
 if(!are_buffers_filled)
   initializeBuffers(viewer);
 attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
 program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
 program->bind();
 vaos[Edges]->bind();
 program->setAttributeValue("colors", QColor(0,0,0));
 glDrawElements(GL_LINES, idx_edge_data_.size(),
                GL_UNSIGNED_INT, idx_edge_data_.data());
 vaos[Edges]->release();
 program->release();
}

void Scene_surface_mesh_item::draw_points(CGAL::Three::Viewer_interface *) const
{

}

bool
Scene_surface_mesh_item::supportsRenderingMode(RenderingMode m) const
{ return (m == FlatPlusEdges || m == Wireframe || m == Flat || m == Gouraud); }

CGAL::Three::Scene_item::Bbox Scene_surface_mesh_item::bbox() const
{
  SMesh::Property_map<vertex_descriptor, Point> pprop = smesh_->points();
  CGAL::Bbox_3 bbox;

  BOOST_FOREACH(vertex_descriptor vd,vertices(*smesh_))
  {
    bbox = bbox + pprop[vd].bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

bool
Scene_surface_mesh_item::isEmpty() const
{

  return num_vertices(*smesh_)==0;
}

QString Scene_surface_mesh_item::toolTip() const
{
  return QObject::tr("<p>Surface_mesh <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of faces: %4</p>")
    .arg(this->name())
    .arg(num_vertices(*smesh_))
    .arg(num_edges(*smesh_))
    .arg(num_faces(*smesh_))
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

#include "Scene_surface_mesh_item.moc"
