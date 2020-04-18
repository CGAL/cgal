#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_3_to_face_graph.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <CGAL/Timer.h>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Exact_predicates_exact_constructions_kernel    EK;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef CGAL::Polyhedron_3<EK>                    EK_Polyhedron_3;
typedef K::Segment_3                              Segment_3;
typedef CGAL::Delaunay_triangulation_3<K>         Delaunay;

// define point creator
typedef K::Point_3                                Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;

void load_from_file(const char* path,std::vector<Point_3>& points)
{
  std::ifstream infile (path);
  std::size_t nbpt;
  infile >> nbpt;
  points.reserve(nbpt);
  Point_3 p;
  do
  {
    infile >> p;
    points.push_back(p);
  }
  while (--nbpt>0);
}

int main(int argc,char** argv)
{
  std::vector<Point_3> points;
  if (argc==1){
    CGAL::Random_points_in_sphere_3<Point_3, PointCreator> gen(1.0);
    int nbpt=1000000;
    std::copy_n( gen, nbpt, std::back_inserter(points) );
    std::cout << "Using " << 1000000 << " random points in the unit ball\n";
  }
  else{
    load_from_file(argv[1],points);
    std::cout << "Using a model with " << points.size() << " points.\n";
  }

  Polyhedron_3 poly;

  // compute convex hull
  CGAL::Timer time;
  time.start();
  CGAL::convex_hull_3(points.begin(), points.end(), poly);
  time.stop();
  std::cout << "Static " << time.time() <<" "<< poly.size_of_vertices() << std::endl;

  poly.clear();

  time.reset();
  time.start();
  Delaunay T(points.begin(), points.end());
  time.stop();
  std::cout << "Delaunay " << time.time() << std::endl;
  time.start();
  CGAL::convex_hull_3_to_face_graph(T,poly);
  time.stop();
  std::cout << "Delaunay+to_poly " << time.time() <<" "<< poly.size_of_vertices() << std::endl;
  poly.clear();

  time.reset();
  time.start();
  CGAL::convex_hull_incremental_3( points.begin(), points.end(), poly, false);
  time.stop();
  std::cout << "incremental EPIC " << time.time() <<" "<< poly.size_of_vertices() << std::endl;

  EK_Polyhedron_3 poly2;
  std::vector<EK::Point_3> ek_points;
  ek_points.reserve(points.size());
  CGAL::Cartesian_converter<K,EK> convert;
  for (std::vector<K::Point_3>::iterator it=points.begin();it!=points.end();++it){
    ek_points.push_back(convert(*it));
  }
  time.reset();
  time.start();
  CGAL::convex_hull_incremental_3( ek_points.begin(), ek_points.end(), poly2, false);
  time.stop();
  std::cout << "incremental EPEC " << time.time() <<" "<< poly2.size_of_vertices() << std::endl;

  return 0;
}
