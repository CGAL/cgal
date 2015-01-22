#include <iostream>
#include <fstream>
#include <vector>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Timer.h>
#include <boost/foreach.hpp>

// this can produce false negatives
// mesh should be pure triangle
template<typename PolygonMesh, typename Point>
void generate_near_boundary(const PolygonMesh& mesh,
                            std::vector<Point>& points,
                            std::vector<bool>& on_boundary)
{
  CGAL_assertion(CGAL::is_pure_triangle(mesh));

  typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type Ppmap;

  std::size_t exp_size = num_vertices(mesh) + num_faces(mesh) + num_edges(mesh);
  points.reserve(exp_size);
  on_boundary.reserve(exp_size);

  Ppmap ppmap = get(boost::vertex_point, mesh);
  // put vertices
  vertex_descriptor vd;

  BOOST_FOREACH(vertex_descriptor vb, vertices(mesh)) {
    points.push_back(ppmap[vb]);
    on_boundary.push_back(true);
  }
  // sample middle of edges
  BOOST_FOREACH(edge_descriptor eb, edges(mesh)) {
    halfedge_descriptor hd = halfedge(eb, mesh);
    const Point& p0 = ppmap[target(hd, mesh)];
    const Point& p1 = ppmap[target(opposite(hd, mesh), mesh)];
    const Point& m = CGAL::ORIGIN + (((p0 + (p1 - CGAL::ORIGIN)) - CGAL::ORIGIN) / 2.0);

    bool has_on = false;
    has_on |= CGAL::Triangle_3<K>(ppmap[target(hd, mesh)],
                                  ppmap[target(next(hd, mesh), mesh)],
                                  ppmap[target(prev(hd, mesh), mesh)]).has_on(m);
    has_on |= CGAL::Triangle_3<K>(ppmap[target(opposite(hd, mesh), mesh)],
                                  ppmap[target(next(opposite(hd, mesh), mesh), mesh)],
                                  ppmap[target(prev(opposite(hd, mesh), mesh), mesh)]).has_on(m);

    points.push_back(m);
    on_boundary.push_back(has_on);
  }
  // sample middle of facets
  BOOST_FOREACH(face_descriptor fb, faces(mesh)) {
    const Point& p0 = ppmap[target(halfedge(fb, mesh), mesh)];
    const Point& p1 = ppmap[target(next(halfedge(fb, mesh), mesh), mesh)];
    const Point& p2 = ppmap[target(prev(halfedge(fb, mesh), mesh), mesh)];

    const Point& m = CGAL::centroid(p0, p1, p2);
    bool has_on = CGAL::Triangle_3<K>(p0, p1, p2).has_on(m);

    points.push_back(m);
    on_boundary.push_back(has_on);
  }
}

template<class OutputIterator, typename PolygonMesh>
void random_points(const PolygonMesh& mesh,
                   int n,
                   OutputIterator out)
{
  CGAL::Bbox_3 bb = bbox(mesh);
  CGAL::Random rg(1340818006); // seed some value for make it easy to debug

  double grid_dx = bb.xmax() - bb.xmin();
  double grid_dy = bb.ymax() - bb.ymin();
  double grid_dz = bb.zmax() - bb.zmin();

  for (int i = 0; i < n; i++){
    *out++ = Point(bb.xmin() + rg.get_double()* grid_dx,
      bb.ymin() + rg.get_double()* grid_dy,
      bb.zmin() + rg.get_double()* grid_dz);
  }
}

template<typename PolygonMesh, typename Point>
void test(
  const PolygonMesh& mesh,
  const std::vector<Point>& points,
  const std::vector<bool>& on_boundary = std::vector<bool>())
{
  std::cerr << "|V| = " << num_vertices(mesh) << std::endl;

  CGAL::Timer timer;
  timer.start();
  CGAL::Point_inside_polygon_mesh<PolygonMesh, K> inside(mesh);
  std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
  timer.reset();

  int nb_inside = 0;
  int nb_boundary = 0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    CGAL::Bounded_side res = inside(points[i]);

    if (!on_boundary.empty()) {
      CGAL_assertion(on_boundary[i] == (res == CGAL::ON_BOUNDARY));
    }

    if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
  }
  timer.stop();
  std::cerr << "Total query size: " << points.size() << std::endl;
  std::cerr << "  " << nb_inside << " points inside " << std::endl;
  std::cerr << "  " << nb_boundary << " points boundary " << std::endl;
  std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;
  std::cerr << " Queries took " << timer.time() << " sec." << std::endl;
}

template<typename PolygonMesh>
CGAL::Bbox_3 bbox(const PolygonMesh& mesh)
{
  typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type
      ppmap = get(boost::vertex_point, mesh);

  CGAL::Bbox_3 bbox(ppmap[*vertices(mesh).first].bbox());
  BOOST_FOREACH(vertex_descriptor vb, vertices(mesh))
  {
    bbox = bbox + ppmap[vb].bbox();
  }
  return bbox;
}

template<typename Point>
void to_file(const char* file_name, const std::vector<Point>& points)
{
  std::ofstream out(file_name);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();
}

template<typename Point>
void from_file(const char* file_name, std::vector<Point>& points)
{
  std::ifstream in(file_name);
  std::copy(std::istream_iterator<Point>(in),
            std::istream_iterator<Point>(),
            std::back_inserter(points));
  in.close();
}
