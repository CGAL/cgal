#include <vector>
#include <queue>

#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <QObject>
#include <QApplication>
#include <QtDebug>

#include <set>
#include <stack>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/version.h> 

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "triangulate_primitive.h"
#include <CGAL/array.h>

#include <map>
struct Scene_polygon_soup_item_priv{

  typedef Polygon_soup::Polygons::const_iterator Polygons_iterator;
  typedef Kernel::Point_3 Point_3;

  Scene_polygon_soup_item_priv(Scene_polygon_soup_item* parent)
    : soup(0),
      oriented(false)
  {
    item = parent;
    nb_polys = 0;
    nb_lines = 0;
    nb_nm_edges = 0;
  }
  ~Scene_polygon_soup_item_priv()
  {
    if(soup)
    {
      delete soup;
      soup = NULL;
    }
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices(void) const;
  void triangulate_polygon(Polygons_iterator, int ) const;
  mutable QOpenGLShaderProgram *program;


  enum VAOs {
      Flat_Facets=0,
      Smooth_Facets,
      Edges,
      NM_Edges,
      NbOfVaos
  };
  enum VBOs {
      Facets_vertices = 0,
      Facets_normals,
      Edges_vertices,
      NM_Edges_vertices,
      F_Colors,
      V_Colors,
      NbOfVbos
  };
  Polygon_soup* soup;
  bool oriented;
  mutable std::vector<float> positions_poly;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> f_colors;
  mutable std::vector<float> v_colors;
  mutable std::vector<float> normals;
  mutable std::vector<float> positions_nm_lines;
  mutable std::size_t nb_nm_edges;
  mutable std::size_t nb_polys;
  mutable std::size_t nb_lines;
  Scene_polygon_soup_item* item;

};


struct Polyhedron_to_polygon_soup_writer {
  typedef Kernel::Point_3 Point_3;

  Polygon_soup* soup;
  Polygon_soup::Polygon_3 polygon;

  Polyhedron_to_polygon_soup_writer(Polygon_soup* soup) : soup(soup), polygon() {
  }

  void write_header( std::ostream&,
                     std::size_t /* vertices */,
                     std::size_t /* halfedges */,
                     std::size_t /* facets */,
                     bool /* normals */ = false ) {
    soup->clear();
  }

  void write_footer() {
  }

  void write_vertex( const double& x, const double& y, const double& z) {
    soup->points.push_back(Point_3(x, y, z));
  }

  void write_normal( const double& /* x */, const double& /* y */, const double& /* z */) {
  }

  void write_facet_header() {
  }

  void write_facet_begin( std::size_t no) {
    polygon.clear();
    polygon.reserve(no);
  }
  void write_facet_vertex_index( std::size_t index) {
    polygon.push_back(index);
  }
  void write_facet_end() {
    soup->polygons.push_back(polygon);
    polygon.clear();
  }
}; // end struct Polyhedron_to_soup_writer

void
Scene_polygon_soup_item_priv::initializeBuffers(CGAL::Three::Viewer_interface* viewer) const
{
    //vao containing the data for the facets
    {
        program = item->getShaderProgram(Scene_polygon_soup_item::PROGRAM_WITH_LIGHT, viewer);
        program->bind();
        item->vaos[Flat_Facets]->bind();
        item->buffers[Facets_vertices].bind();
        item->buffers[Facets_vertices].allocate(positions_poly.data(),
                            static_cast<int>(positions_poly.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Facets_vertices].release();



        item->buffers[Facets_normals].bind();
        item->buffers[Facets_normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);

        item->buffers[Facets_normals].release();
        if(!f_colors.empty())
        {
          item->buffers[F_Colors].bind();
          item->buffers[F_Colors].allocate(f_colors.data(),
                                     static_cast<int>(f_colors.size()*sizeof(float)));
          program->enableAttributeArray("colors");
          program->setAttributeBuffer("colors",GL_FLOAT,0,3);
          item->buffers[F_Colors].release();
        }
        item->vaos[Flat_Facets]->release();



        item->vaos[Smooth_Facets]->bind();
        item->buffers[Facets_vertices].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Facets_vertices].release();

        item->buffers[Facets_normals].release();
        item->buffers[Facets_normals].bind();
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        item->buffers[Facets_normals].release();
        if(!v_colors.empty())
        {
          item->buffers[V_Colors].bind();
          item->buffers[V_Colors].allocate(v_colors.data(),
                                     static_cast<int>(v_colors.size()*sizeof(float)));
          program->enableAttributeArray("colors");
          program->setAttributeBuffer("colors",GL_FLOAT,0,3);
          item->buffers[V_Colors].release();
        }
        item->vaos[Smooth_Facets]->release();
        program->release();
        nb_polys = positions_poly.size();
        positions_poly.resize(0);
        std::vector<float>(positions_poly).swap(positions_poly);

        normals.resize(0);
        std::vector<float>(normals).swap(normals);
        v_colors.resize(0);
        std::vector<float>(v_colors).swap(v_colors);

    }
    //vao containing the data for the edges
    {
        program = item->getShaderProgram(Scene_polygon_soup_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        item->vaos[Edges]->bind();

        item->buffers[Edges_vertices].bind();
        item->buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Edges_vertices].release();
        program->release();
        item->vaos[Edges]->release();

        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<float>(positions_lines).swap(positions_lines);

    }
    //vao containing the data for the non manifold edges
    {
        program = item->getShaderProgram(Scene_polygon_soup_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        item->vaos[NM_Edges]->bind();
        item->buffers[NM_Edges_vertices].bind();
        item->buffers[NM_Edges_vertices].allocate(positions_nm_lines.data(),
                            static_cast<int>(positions_nm_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[NM_Edges_vertices].release();
        item->vaos[NM_Edges]->release();
        nb_nm_edges = positions_nm_lines.size();
        positions_nm_lines.resize(0);
        std::vector<float> (positions_nm_lines).swap(positions_nm_lines);
    }
    item->are_buffers_filled = true;
}

typedef Polyhedron::Traits Traits;
typedef Polygon_soup::Polygon_3 Facet;

void
Scene_polygon_soup_item_priv::triangulate_polygon(Polygons_iterator pit, int polygon_id) const
{
    //Computes the normal of the facet
    const Point_3& pa = soup->points[pit->at(0)];
    const Point_3& pb = soup->points[pit->at(1)];
    const Point_3& pc = soup->points[pit->at(2)];
    Traits::Vector_3 normal = CGAL::cross_product(pb-pa, pc -pa);
    normal = normal / std::sqrt(normal * normal);
    typedef FacetTriangulator<Polyhedron, Kernel, std::size_t> FT;

    double diagonal;
    if(item->diagonalBbox() != std::numeric_limits<double>::infinity())
      diagonal = item->diagonalBbox();
    else
      diagonal = 0.0;
    std::size_t it = 0;
    std::size_t it_end =pit->size();
    std::vector<FT::PointAndId> pointIds;
    do {
      FT::PointAndId pointId;

      pointId.point = soup->points[pit->at(it)];
      pointId.id = pit->at(it);
      pointIds.push_back(pointId);
    } while( ++it != it_end );
    //detect degenerated faces
    std::vector<FT::PointAndId> pid_stack = pointIds;
    for(std::size_t i = 0; i< pointIds.size(); ++i)
    {
     FT::PointAndId pid = pid_stack.back();
     pid_stack.pop_back();
     BOOST_FOREACH(FT::PointAndId poai, pid_stack)
     {
      if (pid.point== poai.point)
      {
        return;
      }
     }
    }
    FT triangulation(pointIds,normal,diagonal);
    //iterates on the internal faces to add the vertices to the positions
    //and the normals to the appropriate vectors
    for(FT::CDT::Finite_faces_iterator
        ffit = triangulation.cdt->finite_faces_begin(),
        end = triangulation.cdt->finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;

        positions_poly.push_back(ffit->vertex(0)->point().x());
        positions_poly.push_back(ffit->vertex(0)->point().y());
        positions_poly.push_back(ffit->vertex(0)->point().z());
        positions_poly.push_back(1.0);


        positions_poly.push_back(ffit->vertex(1)->point().x());
        positions_poly.push_back(ffit->vertex(1)->point().y());
        positions_poly.push_back(ffit->vertex(1)->point().z());
        positions_poly.push_back(1.0);

        positions_poly.push_back(ffit->vertex(2)->point().x());
        positions_poly.push_back(ffit->vertex(2)->point().y());
        positions_poly.push_back(ffit->vertex(2)->point().z());
        positions_poly.push_back(1.0);


        const Point_3& pa = soup->points[pit->at(0)];
        const Point_3& pb = soup->points[pit->at(1)];
        const Point_3& pc = soup->points[pit->at(2)];

        Kernel::Vector_3 n = CGAL::cross_product(pb-pa, pc -pa);
        n = n / std::sqrt(n * n);
        CGAL::Color color;
        if(!soup->fcolors.empty())
          color = soup->fcolors[polygon_id];
        for(int i=0; i<3; i++)
        {
          normals.push_back(n.x());
          normals.push_back(n.y());
          normals.push_back(n.z());
          if(!soup->fcolors.empty())
          {
            f_colors.push_back((float)color.red()/255);
            f_colors.push_back((float)color.green()/255);
            f_colors.push_back((float)color.blue()/255);
          }
          if(!soup->vcolors.empty())
          {
            CGAL::Color vcolor = soup->vcolors[triangulation.v2v[ffit->vertex(i)]];
            v_colors.push_back((float)vcolor.red()/255);
            v_colors.push_back((float)vcolor.green()/255);
            v_colors.push_back((float)vcolor.blue()/255);
          }
        }
    }
}
void
Scene_polygon_soup_item_priv::compute_normals_and_vertices() const{

    QApplication::setOverrideCursor(Qt::WaitCursor);
    //get the vertices and normals
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();

    typedef Polygon_soup::Polygons::size_type size_type;
    positions_poly.resize(0);
    positions_lines.resize(0);
    normals.resize(0);
    positions_nm_lines.resize(0);
    f_colors.resize(0);
    v_colors.resize(0);
    soup->fill_edges();
    int nb = 0;
    for(Polygons_iterator it = soup->polygons.begin();
        it != soup->polygons.end(); ++it)
    {
        if(it->size()!=3)
        {
            triangulate_polygon(it, nb);
        }
        else{
            const Point_3& pa = soup->points[it->at(0)];
            const Point_3& pb = soup->points[it->at(1)];
            const Point_3& pc = soup->points[it->at(2)];

            Kernel::Vector_3 n = CGAL::cross_product(pb-pa, pc -pa);
            n = n / std::sqrt(n * n);

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());


            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());


            for(size_type i = 0; i < it->size(); ++i)
            {
                const Point_3& p = soup->points[it->at(i)];
                positions_poly.push_back(p.x()+offset.x);
                positions_poly.push_back(p.y()+offset.y);
                positions_poly.push_back(p.z()+offset.z);
                positions_poly.push_back(1.0);
                if(!soup->fcolors.empty())
                {
                  const CGAL::Color color = soup->fcolors[nb];
                    f_colors.push_back((float)color.red()/255);
                    f_colors.push_back((float)color.green()/255);
                    f_colors.push_back((float)color.blue()/255);
                }

                if(!soup->vcolors.empty())
                {
                  const CGAL::Color color = soup->vcolors[it->at(i)];
                  v_colors.push_back((float)color.red()/255);
                  v_colors.push_back((float)color.green()/255);
                  v_colors.push_back((float)color.blue()/255);
                }
            }
        }
        nb++;

        //Lines
        for(size_type i = 0; i < it->size(); ++i)
        {

            const Point_3& pa = soup->points[it->at(i)];
            const Point_3& pb = soup->points[it->at((i+1)%it->size())];
            positions_lines.push_back(pa.x()+offset.x);
            positions_lines.push_back(pa.y()+offset.y);
            positions_lines.push_back(pa.z()+offset.z);
            positions_lines.push_back(1.0);

            positions_lines.push_back(pb.x()+offset.x);
            positions_lines.push_back(pb.y()+offset.y);
            positions_lines.push_back(pb.z()+offset.z);
            positions_lines.push_back(1.0);
        }
    }

    //Non manifold edges
    BOOST_FOREACH(const Polygon_soup::Edge& edge,
                    soup->non_manifold_edges)
    {

        const Point_3& a = soup->points[edge[0]];
        const Point_3& b = soup->points[edge[1]];
        positions_nm_lines.push_back(a.x()+offset.x);
        positions_nm_lines.push_back(a.y()+offset.y);
        positions_nm_lines.push_back(a.z()+offset.z);
        positions_nm_lines.push_back(1.0);

        positions_nm_lines.push_back(b.x()+offset.x);
        positions_nm_lines.push_back(b.y()+offset.y);
        positions_nm_lines.push_back(b.z()+offset.z);
        positions_nm_lines.push_back(1.0);
    }
    QApplication::restoreOverrideCursor();
}


Scene_polygon_soup_item::Scene_polygon_soup_item()
    : Scene_item(Scene_polygon_soup_item_priv::NbOfVbos,Scene_polygon_soup_item_priv::NbOfVaos)
{
  d = new Scene_polygon_soup_item_priv(this);
}

Scene_polygon_soup_item::~Scene_polygon_soup_item()
{
  delete d;
}

Scene_polygon_soup_item* 
Scene_polygon_soup_item::clone() const {
  Scene_polygon_soup_item* new_soup = new Scene_polygon_soup_item();
  new_soup->d->soup = d->soup->clone();
  new_soup->d->oriented = d->oriented;
  return new_soup;
}


bool
Scene_polygon_soup_item::load(std::istream& in)
{
  if (!d->soup) d->soup=new Polygon_soup();
  else d->soup->clear();

  bool result = CGAL::read_OFF(in, d->soup->points, d->soup->polygons,
                               d->soup->fcolors, d->soup->vcolors);
  invalidateOpenGLBuffers();
  return result;
}

void Scene_polygon_soup_item::init_polygon_soup(std::size_t nb_pts, std::size_t nb_polygons){
  if(!d->soup)
    d->soup = new Polygon_soup;
  d->soup->clear();
  d->soup->points.reserve(nb_pts);
  d->soup->polygons.reserve(nb_polygons);
  d->soup->vcolors.resize(0);
  d->soup->fcolors.resize(0);
  d->oriented = false;
}



//#include <CGAL/IO/generic_print_polyhedron.h>
#include <iostream>
template<class PolygonMesh>
void polygon_mesh_to_soup(PolygonMesh& mesh, Polygon_soup& soup)
{
  soup.clear();
  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::type VPMap;
  VPMap vpmap = get(boost::vertex_point, mesh);
  std::map<typename boost::graph_traits<PolygonMesh>::vertex_descriptor, int> vim;
  int index=0;
  //fill points
  for(typename boost::graph_traits<PolygonMesh>::vertex_iterator vit =
      vertices(mesh).begin(); vit != vertices(mesh).end(); ++vit)
  {
    soup.points.push_back(get(vpmap, *vit));
    vim.insert(std::make_pair(*vit, index++));
  }
  //fill triangles
  for(typename boost::graph_traits<PolygonMesh>::face_iterator fit =
      faces(mesh).begin(); fit != faces(mesh).end(); ++fit)
  {
    Polygon_soup::Polygon_3 polygon;
    BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor hd,
                  CGAL::halfedges_around_face(halfedge(*fit, mesh), mesh))
    {
      polygon.push_back(vim[target(hd, mesh)]);
    }
    soup.polygons.push_back(polygon);
  }
  soup.fill_edges();

}

void Scene_polygon_soup_item::load(Scene_polyhedron_item* poly_item) {
  if(!poly_item) return;
  if(!poly_item->polyhedron()) return;

  if(!d->soup)
    d->soup = new Polygon_soup;
  polygon_mesh_to_soup(*poly_item->polyhedron(), *d->soup);
  invalidateOpenGLBuffers();
}

void Scene_polygon_soup_item::load(Scene_surface_mesh_item* sm_item) {
  if(!sm_item) return;
  if(!sm_item->face_graph()) return;

  if(!d->soup)
    d->soup = new Polygon_soup;
  polygon_mesh_to_soup(*sm_item->face_graph(), *d->soup);
  invalidateOpenGLBuffers();
}
void
Scene_polygon_soup_item::setDisplayNonManifoldEdges(const bool b)
{

  d->soup->display_non_manifold_edges = b;
}

bool
Scene_polygon_soup_item::displayNonManifoldEdges() const {

  return d->soup->display_non_manifold_edges;
}

void Scene_polygon_soup_item::shuffle_orientations()
{
  for(Polygon_soup::size_type i = 0, end = d->soup->polygons.size();
      i < end; ++i)
  {
    if(std::rand() % 2 == 0) d->soup->inverse_orientation(i);
  }
  invalidateOpenGLBuffers();
}

void Scene_polygon_soup_item::inside_out()
{
  for(Polygon_soup::size_type i = 0, end = d->soup->polygons.size();
      i < end; ++i)
  {
    d->soup->inverse_orientation(i);
  }
  invalidateOpenGLBuffers();
}

bool 
Scene_polygon_soup_item::orient()
{

  if(isEmpty() || d->oriented)
    return true; // nothing to do
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->oriented=true;

  //first skip degenerate polygons
  Polygon_soup::Polygons valid_polygons;
  valid_polygons.reserve(d->soup->polygons.size());
  BOOST_FOREACH(Polygon_soup::Polygon_3& polygon, d->soup->polygons)
  {
    std::set<std::size_t> vids;
    bool to_remove=false;
    BOOST_FOREACH(std::size_t id, polygon)
    {
      if (!vids.insert(id).second){
        to_remove=true;
        break;
      }
    }
    if (!to_remove) valid_polygons.push_back(polygon);
  }
  QApplication::restoreOverrideCursor();
  if (valid_polygons.size()!=d->soup->polygons.size())
    d->soup->polygons.swap(valid_polygons);

  bool res;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  res =  CGAL::Polygon_mesh_processing::
    orient_polygon_soup(d->soup->points, d->soup->polygons);
  QApplication::restoreOverrideCursor();
  return res;
}


bool 
Scene_polygon_soup_item::save(std::ostream& out) const
{

  typedef Polygon_soup::size_type size_type;
  CGAL::File_writer_OFF writer;
  writer.write_header(out,
                      d->soup->points.size(),
                      0,
                      d->soup->polygons.size());
  for(size_type i = 0, end = d->soup->points.size();
      i < end; ++i)
  {
    const Point_3& p = d->soup->points[i];
    writer.write_vertex( p.x(), p.y(), p.z() );
  }
  writer.write_facet_header();
  for(size_type i = 0, end = d->soup->polygons.size();
      i < end; ++i)
  {
    const Polygon_soup::Polygon_3& polygon = d->soup->polygons[i];
    const size_type size = polygon.size();
    writer.write_facet_begin(size);
    for(size_type j = 0; j < size; ++j) {
      writer.write_facet_vertex_index(polygon[j]);
    }
    writer.write_facet_end();
  }
  writer.write_footer();

  return (bool) out;
}

bool 
Scene_polygon_soup_item::exportAsPolyhedron(Polyhedron* out_polyhedron)
{
  if (!orient())
    return false;

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<Polyhedron>(
    d->soup->points, d->soup->polygons, *out_polyhedron);
  std::size_t rv = CGAL::Polygon_mesh_processing::remove_isolated_vertices(*out_polyhedron);
  if(rv > 0)
    std::cerr << "Ignore isolated vertices: " << rv << std::endl;
  if(out_polyhedron->size_of_vertices() > 0) {
    // Also check whether the consistent orientation is fine
    if(out_polyhedron->is_closed() &&
       !CGAL::Polygon_mesh_processing::is_outward_oriented(*out_polyhedron)) {
      out_polyhedron->inside_out();
    }
    return true;
  }
  return false;
}

bool
Scene_polygon_soup_item::exportAsSurfaceMesh(SMesh *out_surface_mesh)
{
  orient();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh< CGAL::Surface_mesh<Point_3> >(
    d->soup->points, d->soup->polygons, *out_surface_mesh);
  std::size_t rv = CGAL::Polygon_mesh_processing::remove_isolated_vertices(*out_surface_mesh);
  if(rv > 0)
    std::cerr << "Ignore isolated vertices: " << rv << std::endl;
  if(out_surface_mesh->vertices().size() > 0) {
    return true;
  }
  return false;
}
QString 
Scene_polygon_soup_item::toolTip() const
{

  if(!d->soup)
    return QString();

  QString str = QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                     "<i>Polygons soup</i></p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of polygons: %3</p>")
    .arg(this->name())
    .arg(d->soup->points.size())
    .arg(d->soup->polygons.size())
    .arg(this->renderingModeName())
    .arg(this->color().name());
    str += QString("<br />Number of isolated vertices: %1<br />").arg(getNbIsolatedvertices());
    return str;
}

void
Scene_polygon_soup_item::draw(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
     d->compute_normals_and_vertices();
     d->initializeBuffers(viewer);
    }
    if(d->soup == 0) return;
    attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    d->program = getShaderProgram(PROGRAM_WITH_LIGHT);
    d->program->bind();

    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
    {
      vaos[Scene_polygon_soup_item_priv::Flat_Facets]->bind();
      if(d->soup->fcolors.empty())
      {
        d->program->setAttributeValue("colors", this->color());
      }
    }
    else if(renderingMode() == Gouraud)
    {
     vaos[Scene_polygon_soup_item_priv::Smooth_Facets]->bind();
      if(d->soup->vcolors.empty())
        d->program->setAttributeValue("colors", this->color());
    }
    //draw the polygons
    // the third argument is the number of vec4 that will be entered
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->nb_polys/4));
    // Clean-up
    d->program->release();
    vaos[Scene_polygon_soup_item_priv::Smooth_Facets]->release();
    vaos[Scene_polygon_soup_item_priv::Flat_Facets]->release();
  }

void
Scene_polygon_soup_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
      d->compute_normals_and_vertices();
      d->initializeBuffers(viewer);
    }
    if(d->soup == 0) return;
    vaos[Scene_polygon_soup_item_priv::Edges]->bind();
    attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    QColor color = this->color();
    d->program->setAttributeValue("colors", color);
    //draw the points
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(d->nb_lines/4));
    // Clean-up
    d->program->release();
    vaos[Scene_polygon_soup_item_priv::Edges]->release();
}

void
Scene_polygon_soup_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
  {
     d->compute_normals_and_vertices();
     d->initializeBuffers(viewer);
  }
    if(d->soup == 0) return;
    vaos[Scene_polygon_soup_item_priv::Edges]->bind();
    attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    d->program->setAttributeValue("colors", QColor(Qt::black));
    //draw the edges
    viewer->glDrawArrays(GL_LINES, 0,static_cast<GLsizei>( d->nb_lines/4));
    // Clean-up
    d->program->release();
    vaos[Scene_polygon_soup_item_priv::Edges]->release();
    if(displayNonManifoldEdges())
    {
        vaos[Scene_polygon_soup_item_priv::NM_Edges]->bind();
        attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
        d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        d->program->bind();
        QColor c = QColor(Qt::red);
        d->program->setAttributeValue("colors", c);
        //draw the edges
        viewer->glDrawArrays(GL_LINES, 0,static_cast<GLsizei>( d->nb_nm_edges/4));
        // Clean-up
        d->program->release();
        vaos[Scene_polygon_soup_item_priv::NM_Edges]->release();
    }

}

bool
Scene_polygon_soup_item::isEmpty() const {

  return (d->soup == 0 || d->soup->points.empty());
}
void
Scene_polygon_soup_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}

void Scene_polygon_soup_item::compute_bbox() const {

  if (isEmpty())
    return;
  const Point_3& p = *(d->soup->points.begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polygon_soup::Points::const_iterator it = d->soup->points.begin();
      it != d->soup->points.end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void 
Scene_polygon_soup_item::new_vertex(const double& x,
                                    const double& y,
                                    const double& z)
{

    d->soup->points.push_back(Point_3(x, y, z));
}
                               
void 
Scene_polygon_soup_item::new_triangle(const std::size_t i,
                                      const std::size_t j,
                                      const std::size_t k)
{

  Polygon_soup::Polygon_3 new_polygon(3);
  new_polygon[0] = i;
  new_polygon[1] = j;
  new_polygon[2] = k;
  d->soup->polygons.push_back(new_polygon);
}

template <class Point, class Polygon>
void Scene_polygon_soup_item::load(const std::vector<Point>& points, const std::vector<Polygon>& polygons)
{
    if(!d->soup)
        d->soup = new Polygon_soup;
    d->soup->clear();

    /// add points
    d->soup->points.reserve(points.size());
    BOOST_FOREACH(const Point& p, points)
            d->soup->points.push_back( Point_3(p[0], p[1], p[2]) );

    /// add polygons
    std::size_t nb_polygons=polygons.size();
    d->soup->polygons.resize(nb_polygons);
    for(std::size_t i=0; i<nb_polygons; ++i)
        d->soup->polygons[i].assign(polygons[i].begin(), polygons[i].end());

    /// fill non-manifold edges container
    //soup->fill_edges();
    d->oriented = false;
    invalidateOpenGLBuffers();
}
// Force the instanciation of the template function for the types used in the STL_io_plugin. This is needed
// because the d-pointer forbid the definition in the .h for this function.
template SCENE_POLYGON_SOUP_ITEM_EXPORT void Scene_polygon_soup_item::load<CGAL::cpp11::array<double, 3>, CGAL::cpp11::array<int, 3> >
(const std::vector<CGAL::cpp11::array<double, 3> >& points, const std::vector<CGAL::cpp11::array<int, 3> >& polygons);
template SCENE_POLYGON_SOUP_ITEM_EXPORT void Scene_polygon_soup_item::load<CGAL::Epick::Point_3, std::vector<std::size_t> >
(const std::vector<CGAL::Epick::Point_3>& points, const std::vector<std::vector<std::size_t> >& polygons);

// Local Variables:
// c-basic-offset: 4
// End:

const Scene_polygon_soup_item::Points& Scene_polygon_soup_item::points() const { return d->soup->points; }
bool Scene_polygon_soup_item::isDataColored() { return d->soup->fcolors.size()>0 || d->soup->vcolors.size()>0;}
std::vector<CGAL::Color> Scene_polygon_soup_item::getVColors() const{return d->soup->vcolors;}
std::vector<CGAL::Color> Scene_polygon_soup_item::getFColors() const{return d->soup->fcolors;}

void Scene_polygon_soup_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);
  if(d && item == this)
  {
    if(d->soup)
    {
      delete d->soup;
      d->soup=NULL;
    }
  }
}
