#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Linear_k;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>      Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>   Circular_k;
typedef Circular_k::Point_2                             Point_2;
typedef Circular_k::Circle_2                            Circle_2;
typedef Circular_k::Circular_arc_2                      Circular_arc_2;
typedef Circular_k::Line_arc_2                          Line_arc_2;

int main()
{
  std::cout << "What is the probability that two arcs formed by" << std::endl;
  std::cout << "three random counterclockwise-oriented points on" << std::endl;
  std::cout << "an unity square intersect?" << std::endl;
  CGAL::Random_points_in_square_2<Point_2> g(1.0);
  double prob = 0.0;
  for (int i = 0; i < 1000; i++) {
    Point_2 p1, p2, p3, p4, p5, p6;
    p1 = *g++; p2 = *g++; p3 = *g++;
    p4 = *g++; p5 = *g++; p6 = *g++;
    if(Linear_k().orientation_2_object()(Linear_k::Point_2(p1.x(), p1.y()),
                                         Linear_k::Point_2(p2.x(), p2.y()),
                                         Linear_k::Point_2(p3.x(), p3.y())) !=
      CGAL::COUNTERCLOCKWISE) std::swap(p1, p3);
    Circular_arc_2 arc1 = Circular_arc_2(p1, p2, p3);
    if(Linear_k().orientation_2_object()(Linear_k::Point_2(p4.x(), p4.y()),
                                         Linear_k::Point_2(p5.x(), p5.y()),
                                         Linear_k::Point_2(p6.x(), p6.y())) !=
      CGAL::COUNTERCLOCKWISE) std::swap(p4, p6);
    Circular_arc_2 arc2 = Circular_arc_2(p4, p5, p6);
    std::vector< CGAL::Object > res;
    Circular_k().intersect_2_object()(arc1, arc2, std::back_inserter(res));
    prob += (res.size() != 0) ? 1.0 : 0.0;
  }
  std::cout << "The probability is: " << prob/1000.0 << std::endl;
  return 0;
};
