#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <QObject>
#include <QApplication>
#include <QtDebug>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/version.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE 1
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Polygon_2.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "triangulate_primitive.h"
#include <CGAL/array.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <streambuf>
#include <vector>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Point_container Pc;

struct Scene_polygon_soup_item_priv{

  typedef Polygon_soup::Polygons::const_iterator Polygons_iterator;
  typedef EPICK::Point_3 Point_3;

  Scene_polygon_soup_item_priv(Scene_polygon_soup_item* parent)
    : soup(0),
      oriented(false)
  {
    item = parent;
    nb_polys = 0;
    nb_lines = 0;
    nb_nm_edges = 0;
    invalidate_stats();
  }
  ~Scene_polygon_soup_item_priv()
  {
    if(soup)
    {
      delete soup;
      soup = NULL;
    }
  }
  void compute_normals_and_vertices(void) const;
  void triangulate_polygon(Polygons_iterator, int ) const;
  void invalidate_stats();
  void compute_stats();
  mutable QOpenGLShaderProgram *program;


  enum Face_names {
      Flat_facets=0,
      Smooth_facets,
  };
  enum Edge_names {
    Edges = 0,
    NM_edges
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
  bool is_triangle, is_quad, stats_computed;
  double minl, maxl, meanl, midl, mini, maxi, ave;
  std::size_t nb_null_edges, nb_degen_faces;

  Scene_polygon_soup_item* item;

};

typedef Scene_polygon_soup_item_priv Priv;


typedef EPICK Traits;
typedef Polygon_soup::Polygon_3 Facet;

void
Scene_polygon_soup_item_priv::triangulate_polygon(Polygons_iterator pit, int polygon_id) const
{
  const CGAL::qglviewer::Vec off = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  EPICK::Vector_3 offset(off.x,off.y,off.z);

    //Computes the normal of the facet
    Traits::Vector_3 normal = CGAL::NULL_VECTOR;

    //Newell's method
    for (std::size_t i = 0; i < pit->size() ; ++ i){
      const Point_3& pa = soup->points[pit->at(i)];
      const Point_3& pb = soup->points[pit->at((i+1)%pit->size())];
      double x = normal.x() + (pa.y()-pb.y())*(pa.z()+pb.z());
      double y = normal.y() + (pa.z()-pb.z())*(pa.x()+pb.x());
      double z = normal.z() + (pa.x()-pb.x())*(pa.y()+pb.y());
      normal = Traits::Vector_3(x,y,z);
    }
    if (normal == CGAL::NULL_VECTOR) // No normal could be computed, return
      return;

    typedef FacetTriangulator<SMesh, EPICK, std::size_t> FT;

    std::size_t it = 0;
    std::size_t it_end =pit->size();
    std::vector<FT::PointAndId> pointIds;
    do {
      FT::PointAndId pointId;

      pointId.point = soup->points[pit->at(it)]+offset;
      pointId.id = pit->at(it);
      pointIds.push_back(pointId);
    } while( ++it != it_end );
    //detect degenerated faces
    std::vector<FT::PointAndId> pid_stack = pointIds;
    for(std::size_t i = 0; i< pointIds.size(); ++i)
    {
     FT::PointAndId pid = pid_stack.back();
     pid_stack.pop_back();
     for(FT::PointAndId poai : pid_stack)
     {
      if (pid.point== poai.point)
      {
        return;
      }
     }
    }
    FT triangulation(pointIds,normal);
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


        positions_poly.push_back(ffit->vertex(1)->point().x());
        positions_poly.push_back(ffit->vertex(1)->point().y());
        positions_poly.push_back(ffit->vertex(1)->point().z());

        positions_poly.push_back(ffit->vertex(2)->point().x());
        positions_poly.push_back(ffit->vertex(2)->point().y());
        positions_poly.push_back(ffit->vertex(2)->point().z());

        CGAL::Color color;
        if(!soup->fcolors.empty())
          color = soup->fcolors[polygon_id];
        for(int i=0; i<3; i++)
        {
          normals.push_back(normal.x());
          normals.push_back(normal.y());
          normals.push_back(normal.z());
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

    //get the vertices and normals
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();

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

            EPICK::Vector_3 n = CGAL::cross_product(pb-pa, pc -pa);
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

            positions_lines.push_back(pb.x()+offset.x);
            positions_lines.push_back(pb.y()+offset.y);
            positions_lines.push_back(pb.z()+offset.z);
        }
    }

    //Non manifold edges
    for(const Polygon_soup::Edge& edge :
                    soup->non_manifold_edges)
    {

        const Point_3& a = soup->points[edge[0]];
        const Point_3& b = soup->points[edge[1]];
        positions_nm_lines.push_back(a.x()+offset.x);
        positions_nm_lines.push_back(a.y()+offset.y);
        positions_nm_lines.push_back(a.z()+offset.z);

        positions_nm_lines.push_back(b.x()+offset.x);
        positions_nm_lines.push_back(b.y()+offset.y);
        positions_nm_lines.push_back(b.z()+offset.z);
    }

}


Scene_polygon_soup_item::Scene_polygon_soup_item()
{
  d = new Scene_polygon_soup_item_priv(this);
  for(int i = 1; i>=0; --i)
  {
    setTriangleContainer(i,
                         new Tc(Vi::PROGRAM_WITH_LIGHT, false));
    setEdgeContainer(i,
                     new Ec(Vi::PROGRAM_NO_SELECTION, false));
  }
  setPointContainer(0,
                    new Pc(Vi::PROGRAM_NO_SELECTION, false));
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

template<class PolygonMesh>
void polygon_mesh_to_soup(PolygonMesh& mesh, Polygon_soup& soup)
{
  soup.clear();
  CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(mesh, soup.points, soup.polygons);
  soup.fill_edges();
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
  for(Polygon_soup::Polygon_3& polygon : d->soup->polygons)
  {
    std::set<std::size_t> vids;
    bool to_remove=false;
    for(std::size_t id : polygon)
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
Scene_polygon_soup_item::exportAsSurfaceMesh(SMesh *out_surface_mesh)
{
  orient();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh< CGAL::Surface_mesh<Point_3> >(
    d->soup->points, d->soup->polygons, *out_surface_mesh);
  std::size_t rv = CGAL::Polygon_mesh_processing::remove_isolated_vertices(*out_surface_mesh);
  if(rv > 0){
    std::cerr << "Ignore isolated vertices: " << rv << std::endl;
    out_surface_mesh->collect_garbage();
  }
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
                     "<i>Polygon soup</i></p>"
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
    if(d->soup == 0) return;
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

    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
    {

      if(d->soup->fcolors.empty())
        getTriangleContainer(Priv::Flat_facets)->setColor(this->color());
      getTriangleContainer(Priv::Flat_facets)->setAlpha(alpha());
      getTriangleContainer(Priv::Flat_facets)->draw(viewer, d->soup->fcolors.empty());
    }
    else if(renderingMode() == Gouraud)
    {
      if(d->soup->vcolors.empty())
        getTriangleContainer(Priv::Smooth_facets)->setColor(this->color());
      getTriangleContainer(Priv::Smooth_facets)->setAlpha(alpha());
      getTriangleContainer(Priv::Smooth_facets)->draw(viewer, d->soup->vcolors.empty());
    }
  }

void
Scene_polygon_soup_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const {

    if(d->soup == 0) return;
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
    getPointContainer(0)->setColor(this->color());
    getPointContainer(0)->draw(viewer, true);
}

void
Scene_polygon_soup_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
    if(d->soup == 0) return;
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
    getEdgeContainer(Priv::Edges)->setColor(QColor(Qt::black));
    getEdgeContainer(Priv::Edges)->draw(viewer, true);
    if(displayNonManifoldEdges())
    {
      getEdgeContainer(Priv::NM_edges)->setColor(QColor(Qt::red));
      //draw the edges
      getEdgeContainer(Priv::NM_edges)->draw(viewer, true);
    }
}

bool
Scene_polygon_soup_item::isEmpty() const {

  return (d->soup == 0 || d->soup->points.empty());
}
void
Scene_polygon_soup_item::invalidateOpenGLBuffers()
{
    compute_bbox();
    for(int i=0; i<2; ++i)
    {
      getTriangleContainer(i)->reset_vbos(ALL);
      getEdgeContainer(i)->reset_vbos(ALL);
    }
    getPointContainer(0)->reset_vbos(ALL);
    setBuffersFilled(false);
    d->invalidate_stats();
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
  setBbox(Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax()));
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
    for(const Point& p : points)
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

template <class Point, class Polygon>
void Scene_polygon_soup_item::load(const std::vector<Point>& points, const std::vector<Polygon>& polygons,
                                   const std::vector<CGAL::Color>& fcolors,
                                   const std::vector<CGAL::Color>& vcolors)
{
    load (points, polygons);

    d->soup->fcolors.reserve (fcolors.size());
    std::copy (fcolors.begin(), fcolors.end(), std::back_inserter (d->soup->fcolors));

    d->soup->vcolors.reserve (vcolors.size());
    std::copy (vcolors.begin(), vcolors.end(), std::back_inserter (d->soup->vcolors));
}
// Force the instanciation of the template function for the types used in the STL_io_plugin. This is needed
// because the d-pointer forbid the definition in the .h for this function.
template SCENE_POLYGON_SOUP_ITEM_EXPORT void Scene_polygon_soup_item::load<std::array<double, 3>, std::array<int, 3> >
(const std::vector<std::array<double, 3> >& points, const std::vector<std::array<int, 3> >& polygons);
template SCENE_POLYGON_SOUP_ITEM_EXPORT void Scene_polygon_soup_item::load<CGAL::Epick::Point_3, std::vector<std::size_t> >
(const std::vector<CGAL::Epick::Point_3>& points, const std::vector<std::vector<std::size_t> >& polygons);
template SCENE_POLYGON_SOUP_ITEM_EXPORT void Scene_polygon_soup_item::load<CGAL::Epick::Point_3, std::vector<std::size_t> >
(const std::vector<CGAL::Epick::Point_3>& points, const std::vector<std::vector<std::size_t> >& polygons,
 const std::vector<CGAL::Color>& fcolors,
 const std::vector<CGAL::Color>& vcolors);

// Local Variables:
// c-basic-offset: 4
// End:

const Scene_polygon_soup_item::Points& Scene_polygon_soup_item::points() const { return d->soup->points; }
const Scene_polygon_soup_item::Polygons& Scene_polygon_soup_item::polygons() const { return d->soup->polygons; }
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

const Polygon_soup::Edges&
Scene_polygon_soup_item::non_manifold_edges() const
{
  return d->soup->non_manifold_edges;
}

void Scene_polygon_soup_item::initializeBuffers(Viewer_interface *v) const
{
  getTriangleContainer(Priv::Flat_facets)->initializeBuffers(v);
  getTriangleContainer(Priv::Flat_facets)->setFlatDataSize(d->nb_polys);

  getTriangleContainer(Priv::Smooth_facets)->initializeBuffers(v);
  getTriangleContainer(Priv::Smooth_facets)->setFlatDataSize(d->nb_polys);

  getEdgeContainer(Priv::Edges)->initializeBuffers(v);
  getEdgeContainer(Priv::Edges)->setFlatDataSize(d->nb_lines);

  getEdgeContainer(Priv::NM_edges)->initializeBuffers(v);
  getEdgeContainer(Priv::NM_edges)->setFlatDataSize(d->nb_nm_edges);

  getPointContainer(0)->initializeBuffers(v);
  getPointContainer(0)->setFlatDataSize(d->nb_lines);

  d->normals.resize(0);
  d->positions_poly.resize(0);
  d->normals.shrink_to_fit();
  d->positions_poly.shrink_to_fit();
  d->v_colors.resize(0);
  d->v_colors.shrink_to_fit();
  d->positions_lines.resize(0);
  d->positions_lines.shrink_to_fit();
  d->positions_nm_lines.resize(0);
  d->positions_nm_lines.shrink_to_fit();
}

void Scene_polygon_soup_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->compute_normals_and_vertices();

  getTriangleContainer(Priv::Flat_facets)->allocate(
        Tc::Flat_vertices,
        d->positions_poly.data(),
        static_cast<int>(d->positions_poly.size()*sizeof(float)));

  getTriangleContainer(Priv::Flat_facets)->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));
  if(!d->f_colors.empty())
  {
    getTriangleContainer(Priv::Flat_facets)->allocate(
          Tc::FColors,
          d->f_colors.data(),
          static_cast<int>(d->f_colors.size()*sizeof(float)));
  }
  getTriangleContainer(Priv::Smooth_facets)->allocate(
        Tc::Flat_vertices,
        d->positions_poly.data(),
        static_cast<int>(d->positions_poly.size()*sizeof(float)));

  getTriangleContainer(Priv::Smooth_facets)->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));

  if(!d->v_colors.empty())
  {
    getTriangleContainer(Priv::Smooth_facets)->allocate(
          Tc::VColors,
          d->v_colors.data(),
          static_cast<int>(d->v_colors.size()*sizeof(float)));
  }

  d->nb_polys = d->positions_poly.size();

  getEdgeContainer(Priv::Edges)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));

  getPointContainer(0)->allocate(
        Pc::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));


  getEdgeContainer(Priv::NM_edges)->allocate(
        Ec::Vertices,
        d->positions_nm_lines.data(),
        static_cast<int>(d->positions_nm_lines.size()*sizeof(float)));

  d->nb_nm_edges = d->positions_nm_lines.size();
  d->nb_lines = d->positions_lines.size();

  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}

void Scene_polygon_soup_item::repair(bool erase_dup, bool req_same_orientation)
{
  QApplication::setOverrideCursor(Qt::BusyCursor);
  CGAL::Polygon_mesh_processing::repair_polygon_soup(
        d->soup->points,
        d->soup->polygons,
        CGAL::Polygon_mesh_processing::parameters::
        erase_all_duplicates(erase_dup)
        .require_same_orientation(req_same_orientation));
  QApplication::restoreOverrideCursor();

 // CGAL::Three::Three::information(
}

CGAL::Three::Scene_item::Header_data Scene_polygon_soup_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories

  data.categories.append(std::pair<QString,int>(QString("Vertices"),1));
  data.categories.append(std::pair<QString,int>(QString("Polygons"),4));
  data.categories.append(std::pair<QString,int>(QString("Edges"),6));
  data.categories.append(std::pair<QString,int>(QString("Angles"),3));


  //titles
  data.titles.append(QString("#Points"));

  data.titles.append(QString("#Polygons"));
  data.titles.append(QString("Pure Triangle"));
  data.titles.append(QString("Pure Quad"));
  data.titles.append(QString("#Degenerate Polygons"));

  data.titles.append(QString("#Edges"));
  data.titles.append(QString("Minimum Length"));
  data.titles.append(QString("Maximum Length"));
  data.titles.append(QString("Median Length"));
  data.titles.append(QString("Mean Length"));
  data.titles.append(QString("#Degenerate Edges"));

  data.titles.append(QString("Minimum"));
  data.titles.append(QString("Maximum"));
  data.titles.append(QString("Average"));
  return data;
}

QString Scene_polygon_soup_item::computeStats(int type)
{
  if(!d->stats_computed)
    d->compute_stats();

  switch(type)
  {
  case NB_VERTICES:
    return QString::number(d->soup->points.size());
  case NB_FACETS:
    return QString::number(d->soup->polygons.size());
  case NB_EDGES:
    return QString::number(d->nb_lines/6);

  case NB_DEGENERATED_FACES:
  {
    if(d->is_triangle)
    {
      return QString::number(d->nb_degen_faces);
    }
    else
      return QString("n/a");
  }

  case MIN_LENGTH:
    return QString::number(d->minl);
  case MAX_LENGTH:
    return QString::number(d->maxl);
  case MID_LENGTH:
    return QString::number(d->midl);
  case MEAN_LENGTH:
    return QString::number(d->meanl);
  case NB_NULL_LENGTH:
    return QString::number(d->nb_null_edges);

  case MIN_ANGLE:
    return QString::number(d->mini);
  case MAX_ANGLE:
    return QString::number(d->maxi);
  case MEAN_ANGLE:
    return QString::number(d->ave);

  case IS_PURE_TRIANGLE:
    if(d->is_triangle)
      return QString("yes");
    else
      return QString("no");
  case IS_PURE_QUAD:
    if (d->is_quad)
      return QString("yes");
    else
      return QString("no");
  }
  return QString();
}

void
Scene_polygon_soup_item_priv::
invalidate_stats()
{
  is_triangle = true;
  is_quad = true;
  minl=0;
  maxl=0;
  meanl=0;
  midl=0;
  mini=0;
  maxi=0;
  ave=0;
  nb_null_edges=0;
  nb_degen_faces=0;
  stats_computed = false;
}


void
Scene_polygon_soup_item_priv::compute_stats()
{
  using namespace boost::accumulators;
  accumulator_set< double,
    features< tag::min, tag::max, tag::mean , tag::median> > edges_acc;
  accumulator_set< double,
    features< tag::min, tag::max, tag::mean > > angles_acc;
  double rad_to_deg = 180. / CGAL_PI;


  for(auto poly : soup->polygons)
  {
    if(poly.size() != 3)
      is_triangle = false;
    if(poly.size() != 4)
      is_quad = false;
    for(std::size_t i = 0; i< poly.size(); ++i)
    {
      Polygon_soup::Point_3 a(soup->points[poly[i]]),
          b(soup->points[poly[(i+1)%poly.size()]]),
          c(soup->points[poly[(i+2)%poly.size()]]);
      if (a == b)
        ++nb_null_edges;
      edges_acc(CGAL::sqrt(CGAL::squared_distance(a, b)));
      typename Traits::Vector_3 ba(b, a);
      typename Traits::Vector_3 bc(b, c);
      double cos_angle = (ba * bc)
        / std::sqrt(ba.squared_length() * bc.squared_length());
      if(cos_angle == CGAL_PI || cos_angle == 0)
        ++nb_degen_faces;
      angles_acc(std::acos(cos_angle) * rad_to_deg);
    }
  }

  minl = extract_result< tag::min >(edges_acc);
  maxl = extract_result< tag::max >(edges_acc);
  meanl = extract_result< tag::mean >(edges_acc);
  midl =  extract_result< tag::median >(edges_acc);
  mini = extract_result< tag::min >(angles_acc);
  maxi = extract_result< tag::max >(angles_acc);
  ave = extract_result< tag::mean >(angles_acc);

  stats_computed = true;
}
