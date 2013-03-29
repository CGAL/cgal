#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/Timer.h>

#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_traits_2<K>          Gt;
typedef Gt::Point_2                                         Point;
typedef Gt::Vector_2                                        Vector;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt>       P2DT2;
typedef CGAL::Delaunay_triangulation_2<Gt>                  DT2;

template <class Dt>
class DT2_inserter {
  Dt t;
public:
  template <class Iterator>
  void insert(Iterator begin, Iterator end) {
    t.insert(begin, end);
  }
};


template <class PT, bool large>
class P2DT2_inserter {
  PT t;
public:
  template <class Iterator>
  void insert(Iterator begin, Iterator end) {
    t.insert(begin, end, large);
  }
};

template <class Inserter>
void test_performance(const std::string &name) {
  // Create point sets
  typedef CGAL::Creator_uniform_2<double,Point>  Creator;
  CGAL::Random rnd(7);
  CGAL::Random_points_in_square_2<Point, Creator> in_square(0.5, rnd);

  CGAL::Timer timer;

  for (int n = 1000; n<=1e6; n+=1000) {
    std::vector<Point> pts;
    for (int i=0 ; i<n ; i++) {
      pts.push_back(*in_square++ + Vector(0.5, 0.5));
    }
    
    Inserter inserter;
    
    timer.start();
    inserter.insert(pts.begin(), pts.end());
    timer.stop();
    std::cout << name << "; " << pts.size() << "; " << timer.time() << std::endl;
  }
}

int main() {
  test_performance<DT2_inserter<DT2> >("Euclidean Delaunay");
  test_performance<P2DT2_inserter<P2DT2, false> >("Periodic Delaunay");
  test_performance<P2DT2_inserter<P2DT2, true>  >("Periodic Delaunay, large point set");

  return 0;
}
