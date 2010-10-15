#include "Scene_polygon_soup.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <QtDebug>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

struct Polygon_soup {
  typedef std::vector<Point_3> Points;
  typedef std::vector<std::size_t> Polygon_3;
  typedef std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t> > Edges_map;
  typedef boost::array<std::size_t, 2> Edge;
  typedef std::vector<Polygon_3> Polygons;
  typedef std::vector<Edge> Edges;
  typedef Polygons::size_type size_type;
  Points points;
  Polygons polygons;
  Edges_map edges;
  Edges non_manifold_edges;
  bool display_non_manifold_edges;

  Polygon_soup* clone() const {
    Polygon_soup* result = new Polygon_soup();
    result->points = points;
    result->polygons = polygons;
    result->edges = edges;
    result->non_manifold_edges = non_manifold_edges;
    result->display_non_manifold_edges = display_non_manifold_edges;
    return result;
  }

  void clear() {
    points.clear();
    polygons.clear();
    edges.clear();
    non_manifold_edges.clear();
  }

  void fill_edges() {
    // Fill edges
    edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        edges[std::make_pair(i0, i1)].insert(i);
//         qDebug() << tr("edges[std::make_pair(%1, %2)].insert(%3). Size=%4")
//           .arg(i0).arg(i1).arg(i).arg(edges[std::make_pair(i0, i1)].size());
      }
    }

    // Fill non-manifold edges
    non_manifold_edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        if( (i0 < i1) && 
            (edges[std::make_pair(i0, i1)].size() +
             edges[std::make_pair(i1, i0)].size() > 2) )
        {
          Edge edge;
          edge[0] = i0;
          edge[1] = i1;
          non_manifold_edges.push_back(edge);
        }
      }
    }
  }

  void inverse_orientation(const std::size_t index) {
    std::reverse(polygons[index].begin(), polygons[index].end());
  }
};

Scene_polygon_soup::Scene_polygon_soup()
  : Scene_item_with_display_list(),
    soup(0),
    oriented(false)
{
}

Scene_polygon_soup::~Scene_polygon_soup()
{
  delete soup;
}

Scene_polygon_soup* 
Scene_polygon_soup::clone() const {
  Scene_polygon_soup* new_soup = new Scene_polygon_soup();
  new_soup->soup = soup->clone();
  new_soup->oriented = oriented;
  return new_soup;
}

bool
Scene_polygon_soup::load(std::istream& in)
{
  if(!soup)
    soup = new Polygon_soup;
  CGAL::File_scanner_OFF scanner(in);
  soup->clear();
  soup->points.resize(scanner.size_of_vertices());
  soup->polygons.resize(scanner.size_of_facets());
  for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i) {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    soup->points[i] = Point_3(x, y, z, w);
    scanner.skip_to_next_vertex( i);
  }
  if(!in)
    return false;

  for (std::size_t i = 0; i < scanner.size_of_facets(); ++i) {
    std::size_t no;
    scanner.scan_facet( no, i);
    soup->polygons[i].resize(no);
    for(std::size_t j = 0; j < no; ++j) {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
      {
        soup->polygons[i][j] = id;
      }
      else
        return false;
    }
  }
  soup->fill_edges();
  oriented = false;
  return in;
}

void
Scene_polygon_soup::setDisplayNonManifoldEdges(const bool b)
{
  soup->display_non_manifold_edges = b;
  changed();
}

bool
Scene_polygon_soup::displayNonManifoldEdges() const {
  return soup->display_non_manifold_edges;
}

void Scene_polygon_soup::shuffle_orientations()
{
  for(Polygon_soup::size_type i = 0, end = soup->polygons.size();
      i < end; ++i)
  {
    if(std::rand() % 2 == 0) soup->inverse_orientation(i);
  }
  soup->fill_edges();
  changed();
}

void Scene_polygon_soup::inside_out()
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
Scene_polygon_soup::orient()
{
  typedef Polygon_soup::Polygons::size_type size_type;
  typedef Polygon_soup::Edges_map Edges;

  if(isEmpty() || this->oriented)
    return true; // nothing to do

  Polygon_soup::Polygons& polygons = soup->polygons;
  Polygon_soup::Edges_map& edges = soup->edges;

  std::vector<bool> oriented;
  std::stack<std::size_t> stack;
  using std::make_pair;

  // no polygon is oriented
  oriented.resize(polygons.size());

  size_type polygon_index = 0;
  bool success = true;

  while (polygon_index != polygons.size()) 
  {
    while ( polygon_index != polygons.size() && oriented[polygon_index] ) {
      ++polygon_index;
    }
    if(polygon_index == polygons.size()) break;

//     qDebug() << tr("Seed %1...\n").arg(polygon_index);
    oriented[polygon_index] = true;
    stack.push(polygon_index);
    while(! stack.empty() )
    {
      const size_type to_be_oriented_index = stack.top();
//       qDebug() << tr("polygon #%1").arg(to_be_oriented_index);
      stack.pop();
      const size_type size = polygons[to_be_oriented_index].size();
      for(size_type ih = 0 ; ih < size ; ++ih) {
        size_type ihp1 = ih+1;
        if(ihp1>=size) ihp1 = 0;
        const std::size_t& i1 = polygons[to_be_oriented_index][ih];
        const std::size_t& i2 = polygons[to_be_oriented_index][ihp1];

//         qDebug() << tr("edge %3-%4 (%1,%2)").arg(i1).arg(i2).arg(ih).arg(ihp1);
        // edge (i1,i2)
        Edges::iterator it_same_orient = edges.find(make_pair(i1, i2));
        // edges (i2,i1)
        Edges::iterator it_other_orient = edges.find(make_pair(i2, i1));

        CGAL_assertion(it_same_orient != edges.end());
        if(it_same_orient->second.size() > 1) {
          if((it_other_orient != edges.end() && it_other_orient->second.size() > 0) ||
             it_same_orient->second.size() > 2) {
            // three polygons at the edge
//             qDebug() << "three polygons at the edge";
            success = false; // non-orientable
          }
          {
            // one neighbor polyhedron, opposite orientation
            size_type index = *(it_same_orient->second.begin());
            if(index == to_be_oriented_index)
              index = *(++it_same_orient->second.begin());
            if(oriented[index]) {
//               qDebug() << tr("neighbor polygon #%1 is already oriented, but in opposite orientation").arg(index);
              success = false; // non-orientable
              continue; // next edge
            }
   
            // reverse the orientation
            const size_type size = polygons[index].size();
            for(size_type j = 0; j < size; ++j) {
              const std::size_t& i0 = polygons[index][j];
              const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
              CGAL_assertion_code(const bool r = )
                edges[std::make_pair(i0, i1)].erase(index);
              CGAL_assertion(r);
            }
            soup->inverse_orientation(index);
            for(size_type j = 0; j < size; ++j) {
              const std::size_t& i0 = polygons[index][j];
              const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
              edges[std::make_pair(i0, i1)].insert(index);
            }
//             qDebug() << tr("inverse the orientation of polygon #%1\n").arg(index);
            oriented[index] = true;
            stack.push(index);
          }
        }
        else if(it_other_orient != edges.end() && it_other_orient->second.size() == 1) {
          // one polygon, same orientation
          const size_type index = *(it_other_orient->second.begin());
          if(oriented[index])
            continue;
          oriented[index] = true;
//           qDebug() << tr("keep the orientation of polygon #%1\n").arg(index);
          stack.push(index);
        }
        else {
//           qDebug() << "else" << it_same_orient->second.size() << it_other_orient->second.size();
          success = false; // non-orientable
        }
      } // end for on all edges of one 
    } // end while loop on the polygons of the connected component
  } // end while loop on all non-oriented polygons remaining 
  return success;
}


bool 
Scene_polygon_soup::save(std::ostream& out) const
{
  typedef Polygon_soup::size_type size_type;
  CGAL::File_writer_OFF writer(out);
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

  return out;
}

QString 
Scene_polygon_soup::toolTip() const
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
Scene_polygon_soup::direct_draw() const {
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
    ::glGetDoublev(GL_CURRENT_COLOR, current_color);
    ::glColor3d(1., 0., 0.); // red
    
    for(Polygon_soup::size_type 
          i = 0,
          end = soup->non_manifold_edges.size();
        i < end; ++i) 
    {
      const Polygon_soup::Edge& edge = soup->non_manifold_edges[i];
      const Point_3& a = soup->points[edge[0]];
      const Point_3& b = soup->points[edge[1]];
      ::glBegin(GL_LINES);
      ::glVertex3d(a.x(), a.y(), a.z());
      ::glVertex3d(b.x(), b.y(), b.z());
      ::glEnd();
    }
    ::glColor4dv(current_color);
  }
}

bool
Scene_polygon_soup::isEmpty() const {
  return (soup == 0 || soup->points.empty());
}

Scene_polygon_soup::Bbox
Scene_polygon_soup::bbox() const {
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

#include "Scene_polygon_soup.moc"
