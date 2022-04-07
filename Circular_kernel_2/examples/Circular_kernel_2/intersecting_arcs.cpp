#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>
#include <iostream>

typedef CGAL::Exact_circular_kernel_2             Circular_k;

typedef CGAL::Point_2<Circular_k>                 Point_2;
typedef CGAL::Circle_2<Circular_k>                Circle_2;
typedef CGAL::Circular_arc_2<Circular_k>          Circular_arc_2;

template <typename T>
double prob_2() {
  CGAL::Random_points_in_square_2<Point_2> g(1.0);
  double prob = 0.0;
  for (int i = 0; i < 10000; i++) {

    Point_2 p1, p2, p3, p4, p5, p6;
    p1 = *g++; p2 = *g++; p3 = *g++;
    p4 = *g++; p5 = *g++; p6 = *g++;

    // the pi's are points inherited from the Cartesian kernel Point_2, so,
    // the orientation predicate can be called on them
    if(CGAL::orientation(p1, p2, p3) != CGAL::COUNTERCLOCKWISE) std::swap(p1, p3);
    T o1 = T(p1, p2, p3);
    if(CGAL::orientation(p4, p5, p6) != CGAL::COUNTERCLOCKWISE) std::swap(p4, p6);
    T o2 = T(p4, p5, p6);

    typedef typename CGAL::CK2_Intersection_traits<Circular_k, T, T>::type
      Intersection_result;
    std::vector<Intersection_result> res;
    CGAL::intersection(o1, o2, std::back_inserter(res));

    prob += (res.size() != 0) ? 1.0 : 0.0;
  }
  return prob/10000.0;
}

int main()
{
  std::cout << "What is the probability that two arcs formed by" << std::endl;
  std::cout << "three random counterclockwise-oriented points on" << std::endl;
  std::cout << "an unit square intersect? (wait a second please)" << std::endl;
  std::cout << "The probability is: " << prob_2<Circular_arc_2>() <<
    std::endl << std::endl;

  std::cout << "And what about the probability that two circles formed by"
    << std::endl;
  std::cout << "three random counterclockwise-oriented points on" << std::endl;
  std::cout << "an unit square intersect? (wait a second please)" << std::endl;
  std::cout << "The probability is: " << prob_2<Circle_2>() << std::endl;
  return 0;
}
