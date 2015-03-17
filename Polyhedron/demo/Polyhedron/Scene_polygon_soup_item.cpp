#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <QObject>
#include <QtDebug>

#include <set>
#include <stack>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/version.h> 

#include <CGAL/polygon_soup_to_polygon_mesh.h>
#include <CGAL/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

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

Scene_polygon_soup_item::Scene_polygon_soup_item()
  : Scene_item_with_display_list(),
    soup(0),
    oriented(false)
{
}

Scene_polygon_soup_item::~Scene_polygon_soup_item()
{
  delete soup;
}

Scene_polygon_soup_item* 
Scene_polygon_soup_item::clone() const {
  Scene_polygon_soup_item* new_soup = new Scene_polygon_soup_item();
  new_soup->soup = soup->clone();
  new_soup->oriented = oriented;
  return new_soup;
}

bool
Scene_polygon_soup_item::load(std::istream& in)
{
  if (!soup) soup=new Polygon_soup();
  else soup->clear();
  return CGAL::read_OFF(in, soup->points, soup->polygons);
}

void Scene_polygon_soup_item::init_polygon_soup(std::size_t nb_pts, std::size_t nb_polygons){
  if(!soup)
    soup = new Polygon_soup;
    soup->clear();
  soup->points.reserve(nb_pts);
  soup->polygons.reserve(nb_polygons);
  oriented = false;
}

void Scene_polygon_soup_item::finalize_polygon_soup(){ soup->fill_edges(); }

#include <CGAL/IO/generic_print_polyhedron.h>
#include <iostream>

void Scene_polygon_soup_item::load(Scene_polyhedron_item* poly_item) {
  if(!poly_item) return;
  if(!poly_item->polyhedron()) return;

  if(!soup)
    soup = new Polygon_soup;

  Polyhedron_to_polygon_soup_writer writer(soup);
  CGAL::generic_print_polyhedron(std::cerr,
                                 *poly_item->polyhedron(),
                                 writer);
  emit changed();
}

void
Scene_polygon_soup_item::setDisplayNonManifoldEdges(const bool b)
{
  soup->display_non_manifold_edges = b;
  changed();
}

bool
Scene_polygon_soup_item::displayNonManifoldEdges() const {
  return soup->display_non_manifold_edges;
}

void Scene_polygon_soup_item::shuffle_orientations()
{
  for(Polygon_soup::size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    if(std::rand() % 2 == 0) soup->inverse_orientation(i);
  }
  soup->fill_edges();
  changed();
}

void Scene_polygon_soup_item::inside_out()
{
  for(Polygon_soup::size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    soup->inverse_orientation(i);
  }
  soup->fill_edges();
  changed();
}

bool 
Scene_polygon_soup_item::orient()
{
  if(isEmpty() || oriented)
    return true; // nothing to do
  oriented=true;
  return CGAL::Polygon_mesh_processing::
    orient_polygon_soup(soup->points, soup->polygons);
}


bool 
Scene_polygon_soup_item::save(std::ostream& out) const
{
  typedef Polygon_soup::size_type size_type;
  CGAL::File_writer_OFF writer;
  writer.write_header(out,
                      soup->points.size(),
                      0,
                      soup->polygons.size());
  for(size_type i = 0, end = soup->points.size();
      i < end; ++i)
  {
    const Point_3& p = soup->points[i];
    writer.write_vertex( p.x(), p.y(), p.z() );
  }
  writer.write_facet_header();
  for(size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    const Polygon_soup::Polygon_3& polygon = soup->polygons[i]; 
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
  orient();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<Polyhedron>(
    soup->points, soup->polygons, *out_polyhedron);

  if(out_polyhedron->size_of_vertices() > 0) {
    // Also check whether the consistent orientation is fine
    if(!CGAL::Polygon_mesh_processing::is_outward_oriented(*out_polyhedron)) {
      out_polyhedron->inside_out();
    }
    return true;
  }
  return false;
}

QString 
Scene_polygon_soup_item::toolTip() const
{
  if(!soup)
    return QString();

  return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                     "<i>Polygons soup</i></p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of polygons: %3</p>")
    .arg(this->name())
    .arg(soup->points.size())
    .arg(soup->polygons.size())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

void
Scene_polygon_soup_item::direct_draw() const {
  typedef Polygon_soup::Polygons::const_iterator Polygons_iterator;
  typedef Polygon_soup::Polygons::size_type size_type;
  for(Polygons_iterator it = soup->polygons.begin();
      it != soup->polygons.end(); ++it)
  {
    const Point_3& pa = soup->points[it->at(0)];
    const Point_3& pb = soup->points[it->at(1)];
    const Point_3& pc = soup->points[it->at(2)];

    Kernel::Vector_3 n = CGAL::cross_product(pb-pa, pc -pa);
    n = n / std::sqrt(n * n);
    
    ::glBegin(GL_POLYGON);  
    ::glNormal3d(n.x(),n.y(),n.z());

    for(size_type i = 0; i < it->size(); ++i) {
      const Point_3& p = soup->points[it->at(i)];
      ::glVertex3d(p.x(),p.y(),p.z());
    }
    ::glEnd();
  }
  if(soup->display_non_manifold_edges) {
    double current_color[4];
    GLboolean lightning;
    ::glGetDoublev(GL_CURRENT_COLOR, current_color);
    ::glColor3d(1., 0., 0.); // red
    ::glGetBooleanv(GL_LIGHTING, &lightning);
    ::glDisable(GL_LIGHTING);

    BOOST_FOREACH(const Polygon_soup::Edge& edge,
                  soup->non_manifold_edges)
    {
      const Point_3& a = soup->points[edge[0]];
      const Point_3& b = soup->points[edge[1]];
      ::glBegin(GL_LINES);
      ::glVertex3d(a.x(), a.y(), a.z());
      ::glVertex3d(b.x(), b.y(), b.z());
      ::glEnd();
    }
    if(lightning) glEnable(GL_LIGHTING);
    ::glColor4dv(current_color);
  }
}

void
Scene_polygon_soup_item::draw_points() const {
  if(soup == 0) return;
  ::glBegin(GL_POINTS);
  for(Polygon_soup::Points::const_iterator pit = soup->points.begin(),
        end = soup->points.end();
      pit != end; ++pit)
  {
    ::glVertex3d(pit->x(), pit->y(), pit->z());
  }
  ::glEnd();
}

bool
Scene_polygon_soup_item::isEmpty() const {
  return (soup == 0 || soup->points.empty());
}

Scene_polygon_soup_item::Bbox
Scene_polygon_soup_item::bbox() const {
  const Point_3& p = *(soup->points.begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polygon_soup::Points::const_iterator it = soup->points.begin();
      it != soup->points.end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void 
Scene_polygon_soup_item::new_vertex(const double& x,
                               const double& y,
                               const double& z)
{
  soup->points.push_back(Point_3(x, y, z));
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
  soup->polygons.push_back(new_polygon);
}
                               
#include "Scene_polygon_soup_item.moc"
