#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3                              Segment_3;

// define point creator
typedef K::Point_3                                Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;

//a functor computing the plane containing a triangular facet
struct Plane_from_facet {
  Polyhedron_3::Plane_3 operator()(Polyhedron_3::Facet& f) {
      Polyhedron_3::Halfedge_handle h = f.halfedge();
      return Polyhedron_3::Plane_3( h->vertex()->point(),
                                    h->next()->vertex()->point(),
                                    h->opposite()->vertex()->point());
  }
};


int main()
{
  CGAL::Random_points_in_sphere_3<Point_3, PointCreator> gen(100.0);

  // generate 250 points randomly on a sphere of radius 100.0
  // and copy them to a vector
  std::vector<Point_3> points;
  CGAL::cpp11::copy_n( gen, 250, std::back_inserter(points) );

  // define polyhedron to hold convex hull
  Polyhedron_3 poly;
  
  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
  
  // assign a plane equation to each polyhedron facet using functor Plane_from_facet
  std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),Plane_from_facet());


  return 0;
}
