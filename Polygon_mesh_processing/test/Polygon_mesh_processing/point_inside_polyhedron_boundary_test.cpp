#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_inside_polygon_mesh_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Timer.h>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;


typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::const_type Ppmap;

// this can produce false negatives
// poly should be pure triangle
void generate_near_boundary(const Polyhedron& poly, std::vector<Point>& points, std::vector<bool>& on_boundary) {
  CGAL_assertion(poly.is_pure_triangle());

  std::size_t exp_size = num_vertices(poly) + num_faces(poly) + num_edges(poly);
  points.reserve(exp_size);
  on_boundary.reserve(exp_size);

  Ppmap ppmap = get(boost::vertex_point,poly);
  // put vertices
  vertex_descriptor vd;

  BOOST_FOREACH(vertex_descriptor vb, vertices(poly)) {
    points.push_back(ppmap[vb]);
    on_boundary.push_back(true);
  }
  // sample middle of edges
  BOOST_FOREACH(edge_descriptor eb, edges(poly)) {
    halfedge_descriptor hd = halfedge(eb,poly);
    const Point& p0 = ppmap[target(hd,poly)];
    const Point& p1 = ppmap[target(opposite(hd,poly),poly)];
    const Point& m = CGAL::ORIGIN + (((p0 + (p1 - CGAL::ORIGIN)) - CGAL::ORIGIN ) / 2.0);

    bool has_on = false;
    has_on |= CGAL::Triangle_3<K>(ppmap[target(hd,poly)], ppmap[target(next(hd,poly),poly)],
                                  ppmap[target(prev(hd,poly),poly)]).has_on(m);
    has_on |= CGAL::Triangle_3<K>(ppmap[target(opposite(hd,poly),poly)], ppmap[target(next(opposite(hd,poly),poly),poly)], 
                                  ppmap[target(prev(opposite(hd,poly),poly),poly)]).has_on(m);
    
    points.push_back(m);
    on_boundary.push_back(has_on);
  }
  // sample middle of facets
  BOOST_FOREACH(face_descriptor fb, faces(poly)) {
    const Point& p0 = ppmap[target(halfedge(fb,poly),poly)];
    const Point& p1 = ppmap[target(next(halfedge(fb,poly),poly),poly)];
    const Point& p2 = ppmap[target(prev(halfedge(fb,poly),poly),poly)];

    const Point& m = CGAL::centroid(p0, p1, p2);
    bool has_on = CGAL::Triangle_3<K>(p0, p1, p2).has_on(m);

    points.push_back(m);
    on_boundary.push_back(has_on);
  }
}

template<class OutputIterator>
void random_points(CGAL::Bbox_3 bbox, int n, OutputIterator out)
{
  CGAL::Random rg(1340818006); // seed some value for make it easy to debug

  double grid_dx = bbox.xmax() - bbox.xmin();
  double grid_dy = bbox.ymax() - bbox.ymin();
  double grid_dz = bbox.zmax() - bbox.zmin();

  for(int i=0; i < n; i++){
    *out++ = Point(bbox.xmin() + rg.get_double()* grid_dx, 
                   bbox.ymin() + rg.get_double()* grid_dy,
                   bbox.zmin() + rg.get_double()* grid_dz);
  }
}

void test(
  const Polyhedron& poly,
  const std::vector<Point>& points,
  const std::vector<bool>& on_boundary = std::vector<bool>())
{
  std::cerr << "|V| = " << num_vertices(poly) << std::endl;

  CGAL::Timer timer;
  timer.start();
  CGAL::Point_inside_polygon_mesh<Polyhedron, K> inside(poly);
  std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
  timer.reset();  

  int nb_inside = 0;
  int nb_boundary = 0;
  for(std::size_t i = 0; i < points.size(); ++i) {
    CGAL::Bounded_side res = inside(points[i]);

    if(!on_boundary.empty()) {
      CGAL_assertion(on_boundary[i] == (res == CGAL::ON_BOUNDARY));
    }

    if(res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    if(res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
  }
  timer.stop();
  std::cerr << "Total query size: " << points.size() << std::endl;
  std::cerr << "  " << nb_inside << " points inside " << std::endl;
  std::cerr << "  " << nb_boundary << " points boundary " << std::endl;
  std::cerr << "  " << points.size() - nb_inside -nb_boundary << " points outside "  << std::endl;
  std::cerr << " Queries took " << timer.time() << " sec." << std::endl; 
}

CGAL::Bbox_3 bbox(const Polyhedron& polyhedron) {
  Ppmap ppmap = get(boost::vertex_point,polyhedron);
  CGAL::Bbox_3 bbox(ppmap[*vertices(polyhedron).first].bbox());
  BOOST_FOREACH(vertex_descriptor vb, vertices(polyhedron)) {
    bbox = bbox + ppmap[vb].bbox();
  }
  return bbox;
}

void to_file(const char* file_name, const std::vector<Point>& points) {
  std::ofstream out(file_name);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();
}

void from_file(const char* file_name, std::vector<Point>& points) {
  std::ifstream in(file_name);
  std::copy(std::istream_iterator<Point>(in),
    std::istream_iterator<Point>(),
    std::back_inserter(points));
  in.close();
}

int main(int, char** argv) {
  std::ifstream input(argv[1]);
  Polyhedron poly;
  if ( !input || !(input >> poly) || poly.is_empty() ){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  std::vector<Point> points;

  std::vector<bool> on_boundary;
  generate_near_boundary(poly, points, on_boundary);
  test(poly, points, on_boundary);

  points.clear();
  const int nb_query = (int)1.e6;
  points.reserve(nb_query);
  random_points(bbox(poly), nb_query, back_inserter(points));
  test(poly, points);

  //test compilation of constructor from AABB_tree
  
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> FGTP;
  typedef CGAL::AABB_traits<K, FGTP>    AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;

  AABB_tree tree(faces(poly).first, faces(poly).second, poly);
  CGAL::Point_inside_polygon_mesh<Polyhedron, K> inside_test(tree);

}
