#include "Scene_surface_mesh_item.h"

#include <queue>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::halfedge_descriptor Halfedge_descriptor;
struct Face_info {
  Halfedge_descriptor e[3];
  bool is_external;
};
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Scene_surface_mesh_item::Kernel>                       P_traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_descriptor,P_traits>                                 Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits >                                           Fb1;
typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>                                                Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                                               TDS;
typedef CGAL::Exact_intersections_tag                                                                             Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>                                           CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>                                                           CDT;



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
  floated = false;
  checkFloat();
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->add_property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, Kernel::Vector_3 > fnormals =
      smesh_->add_property_map<face_descriptor, Kernel::Vector_3 >("v:normal").first;
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
        idx_data_.push_back(im[source(hd, *smesh_)]);
      }
    }
    else if(is_quad(halfedge(fd,*smesh_),*smesh_))
    {
      halfedge_descriptor hd = halfedge(fd,*smesh_);
      //1st half
        idx_data_.push_back(im[source(hd, *smesh_)]);
        idx_data_.push_back(im[source(next(hd, *smesh_), *smesh_)]);
        idx_data_.push_back(im[source(next(next(hd, *smesh_), *smesh_), *smesh_)]);

        //2nd half
        idx_data_.push_back(im[source(hd, *smesh_)]);
        idx_data_.push_back(im[source(next(next(hd, *smesh_), *smesh_), *smesh_)]);
        idx_data_.push_back(im[source(prev(hd, *smesh_), *smesh_)]);
    }
    else
    {
      triangulate_facet(fd, &fnormals, 0, &im, true);
    }
  }

  idx_edge_data_.reserve(num_edges(*smesh_) * 2);
  BOOST_FOREACH(edge_descriptor ed, edges(*smesh_))
  {
    idx_edge_data_.push_back(im[source(ed, *smesh_)]);
    idx_edge_data_.push_back(im[target(ed, *smesh_)]);
  }

  has_vcolors = false;
  has_fcolors = false;
  compute_elements();
}

Scene_surface_mesh_item*
Scene_surface_mesh_item::clone() const
{ return new Scene_surface_mesh_item(*this); }

void Scene_surface_mesh_item::addFlatData(Point p, Kernel::Vector_3 n, CGAL::Color *c) const
{

  flat_vertices.push_back((gl_data)p.x());
  flat_vertices.push_back((gl_data)p.y());
  flat_vertices.push_back((gl_data)p.z());

  flat_normals.push_back((gl_data)n.x());
  flat_normals.push_back((gl_data)n.y());
  flat_normals.push_back((gl_data)n.z());

  if(c != NULL)
  {
    f_colors.push_back((float)c->red()/255);
    f_colors.push_back((float)c->green()/255);
    f_colors.push_back((float)c->blue()/255);
  }
}

void Scene_surface_mesh_item::compute_elements()
{
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
        flat_vertices.push_back((gl_data)p.x());
        flat_vertices.push_back((gl_data)p.y());
        flat_vertices.push_back((gl_data)p.z());
        Kernel::Vector_3 n = fnormals[fd];
        flat_normals.push_back((gl_data)n.x());
        flat_normals.push_back((gl_data)n.y());
        flat_normals.push_back((gl_data)n.z());

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

      hd = halfedge(next(halfedge(fd, *smesh_),*smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = halfedge(next(next(halfedge(fd, *smesh_),*smesh_), *smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);
      //2nd half
      hd = halfedge(fd, *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = halfedge(next(next(halfedge(fd, *smesh_),*smesh_), *smesh_), *smesh_);
      addFlatData(positions[source(hd, *smesh_)]
          ,fnormals[fd]
          ,c);

      hd = halfedge(prev(halfedge(fd, *smesh_), *smesh_), *smesh_);
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
      smooth_vertices.push_back((gl_data)p.x());
      smooth_vertices.push_back((gl_data)p.y());
      smooth_vertices.push_back((gl_data)p.z());

      Kernel::Vector_3 n = vnormals[vd];
      smooth_normals.push_back((gl_data)n.x());
      smooth_normals.push_back((gl_data)n.y());
      smooth_normals.push_back((gl_data)n.z());

    }
  }
}
void Scene_surface_mesh_item::initializeBuffers(CGAL::Three::Viewer_interface* viewer)const
{
  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
    smesh_->points();
  SMesh::Property_map<vertex_descriptor, Kernel::Vector_3 > vnormals =
    smesh_->property_map<vertex_descriptor, Kernel::Vector_3 >("v:normal").first;
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
  if(has_fcolors)
  {
    buffers[FColors].bind();
    buffers[FColors].allocate(f_colors.data(),
                             static_cast<int>(f_colors.size()*sizeof(gl_data)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",GL_DATA,0,3);
    buffers[FColors].release();
  }
  vaos[Flat_facets]->release();

  //vao containing the data for the smooth facets
  vaos[Smooth_facets]->bind();
  buffers[Smooth_vertices].bind();
  if(!floated)
  buffers[Smooth_vertices].allocate(positions.data(),
                             static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  else
    buffers[Smooth_vertices].allocate(smooth_vertices.data(),
                               static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_DATA,0,3);
  buffers[Smooth_vertices].release();


  buffers[Smooth_normals].bind();
  if(!floated)
    buffers[Smooth_normals].allocate(vnormals.data(),
                                     static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  else
    buffers[Smooth_normals].allocate(smooth_normals.data(),
                              static_cast<int>(num_vertices(*smesh_)*3*sizeof(gl_data)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals",GL_DATA,0,3);
  buffers[Smooth_normals].release();
  if(has_vcolors)
  {
    buffers[VColors].bind();
    buffers[VColors].allocate(v_colors.data(),
                             static_cast<int>(v_colors.size()*sizeof(gl_data)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",GL_DATA,0,3);
    buffers[VColors].release();
  }
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
  glShadeModel(GL_SMOOTH);
  if(!are_buffers_filled)
    initializeBuffers(viewer);
  attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
  program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
  program->bind();

  if(renderingMode() == Gouraud)
  {
    vaos[Smooth_facets]->bind();
    if(is_selected)
      program->setAttributeValue("is_selected", true);
    else
      program->setAttributeValue("is_selected", false);
      if(!has_vcolors)
        program->setAttributeValue("colors", this->color());
    glDrawElements(GL_TRIANGLES, idx_data_.size(),
                   GL_UNSIGNED_INT, idx_data_.data());
    vaos[Smooth_facets]->release();
  }
  else
  {
    vaos[Flat_facets]->bind();
    program->setAttributeValue("colors", this->color());
    if(is_selected)
      program->setAttributeValue("is_selected", true);
    else
      program->setAttributeValue("is_selected", false);
    if(!has_fcolors)
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
 if(is_selected)
   program->setAttributeValue("is_selected", true);
 else
   program->setAttributeValue("is_selected", false);
 glDrawElements(GL_LINES, idx_edge_data_.size(),
                GL_UNSIGNED_INT, idx_edge_data_.data());
 vaos[Edges]->release();
 program->release();
}

void Scene_surface_mesh_item::draw_points(CGAL::Three::Viewer_interface *) const
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

void Scene_surface_mesh_item::checkFloat()const
{
#if IS_FLOAT == 1
  floated = true;
#endif
}

void
Scene_surface_mesh_item::triangulate_facet(face_descriptor fd,
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
  P_traits cdt_traits(normal);
  CDT cdt(cdt_traits);

  SMesh::Halfedge_around_face_circulator
      he_circ(halfedge(fd,*smesh_), *smesh_),
      he_circ_end(he_circ);


  // Iterates on the vector of facet handles
  boost::container::flat_map<CDT::Vertex_handle, vertex_descriptor> v2v;
  CDT::Vertex_handle previous, first;
  do {
    CDT::Vertex_handle vh = cdt.insert(smesh_->point(source(*he_circ, *smesh_)));
    if(index)
      v2v.insert(std::make_pair(vh, source(*he_circ, *smesh_)));
    if(first == 0) {
      first = vh;
    }
    vh->info() = *he_circ;
    if(previous != 0 && previous != vh) {
      cdt.insert_constraint(previous, vh);
    }
    previous = vh;
  } while( ++he_circ != he_circ_end );
  cdt.insert_constraint(previous, first);
  // sets mark is_external
  for( CDT::All_faces_iterator
       fit2 = cdt.all_faces_begin(),
       end = cdt.all_faces_end();
       fit2 != end; ++fit2)
  {
    fit2->info().is_external = false;
  }
  //check if the facet is external or internal
  std::queue< CDT::Face_handle> face_queue;
  face_queue.push(cdt.infinite_vertex()->face());
  while(! face_queue.empty() ) {
    CDT::Face_handle fh = face_queue.front();
    face_queue.pop();
    if(fh->info().is_external) continue;
    fh->info().is_external = true;
    for(int i = 0; i <3; ++i) {
      if(!cdt.is_constrained(std::make_pair(fh, i)))
      {
        face_queue.push(fh->neighbor(i));
      }
    }
  }
  //iterates on the internal faces
  for( CDT::Finite_faces_iterator
       ffit = cdt.finite_faces_begin(),
       end = cdt.finite_faces_end();
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
      idx_data_.push_back((*im)[v2v[ffit->vertex(0)]]);
      idx_data_.push_back((*im)[v2v[ffit->vertex(1)]]);
      idx_data_.push_back((*im)[v2v[ffit->vertex(2)]]);
    }

  }
}
#include "Scene_surface_mesh_item.moc"
