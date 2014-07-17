#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

// this can produce false negatives
// poly should be pure triangle
void generate_near_boundary(const Polyhedron& poly, std::vector<Point>& points, std::vector<bool>& on_boundary) {
  CGAL_assertion(poly.is_pure_triangle());

  std::size_t exp_size = poly.size_of_vertices() + poly.size_of_facets() + poly.size_of_halfedges() / 2;
  points.reserve(exp_size);
  on_boundary.reserve(exp_size);

  // put vertices
  for(Polyhedron::Vertex_const_iterator vb = poly.vertices_begin(); vb != poly.vertices_end(); ++vb) {
    points.push_back(vb->point());
    on_boundary.push_back(true);
  }
  // sample middle of edges
  for(Polyhedron::Edge_const_iterator eb = poly.edges_begin(); eb != poly.edges_end(); ++eb) {
    const Point& p0 = eb->vertex()->point();
    const Point& p1 = eb->opposite()->vertex()->point();
    const Point& m = CGAL::ORIGIN + (((p0 + (p1 - CGAL::ORIGIN)) - CGAL::ORIGIN ) / 2.0);

    bool has_on = false;
    has_on |= CGAL::Triangle_3<K>(eb->vertex()->point(), eb->next()->vertex()->point(),
      eb->prev()->vertex()->point()).has_on(m);
    has_on |= CGAL::Triangle_3<K>(eb->opposite()->vertex()->point(), eb->opposite()->next()->vertex()->point(), 
      eb->opposite()->prev()->vertex()->point()).has_on(m);
    
    points.push_back(m);
    on_boundary.push_back(has_on);
  }
  // sample middle of facets
  for(Polyhedron::Facet_const_handle fb = poly.facets_begin(); fb != poly.facets_end(); ++fb) {
    const Point& p0 = fb->halfedge()->vertex()->point();
    const Point& p1 = fb->halfedge()->next()->vertex()->point();
    const Point& p2 = fb->halfedge()->prev()->vertex()->point();

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
  std::cerr << "|V| = " << poly.size_of_vertices() << std::endl;

  CGAL::Timer timer;
  timer.start();
  CGAL::Point_inside_polyhedron_3<Polyhedron, K> inside(poly);
  std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
  timer.reset();  

  int nb_inside = 0;
  int nb_boundary = 0;
  for(std::size_t i = 0; i < points.size(); ++i) {
    CGAL::Bounded_side res = inside(points[i]);

    if(!on_boundary.empty()) {
      assert(on_boundary[i] == (res == CGAL::ON_BOUNDARY));
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
  CGAL::Bbox_3 bbox(polyhedron.vertices_begin()->point().bbox());
  for(Polyhedron::Vertex_const_iterator vb = polyhedron.vertices_begin(); vb != polyhedron.vertices_end(); ++vb) {
    bbox = bbox + vb->point().bbox();
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
  if ( !input || !(input >> poly) || poly.empty() ){
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
}