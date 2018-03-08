
#include <CGAL/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                      CGAL::Hyperbolic_octagon_translation>     Traits;
typedef Traits::Point_2                                                         Point;
typedef Traits::Circle_2                                                        Circle;
typedef Traits::Euclidean_line_2                                                Line;
typedef Traits::Construct_intersection_2                                        Construct_intersection_2;

using CGAL::sqrt;
using std::cout;
using std::endl;

int main(int argc, char** argv) {    

    NT F2(2);

    Line ell1(sqrt(F2), sqrt(F2 - sqrt(F2)), - sqrt(F2 + sqrt(F2)));
    Line ell2(sqrt(F2), sqrt(F2), F2*sqrt(F2 + sqrt(F2)));
    NT sx = -sqrt(F2 + sqrt(F2))*(sqrt(F2) + F2*sqrt(F2 - sqrt(F2)))/(sqrt(F2)*sqrt(F2-sqrt(F2)) - F2);
    NT sy = NT(3)*sqrt(F2 + sqrt(F2))*sqrt(F2)/(sqrt(F2)*sqrt(F2 - sqrt(F2)) - F2);
    Point sol(sx, sy);

    Point p = Construct_intersection_2()(ell1, ell2);
    cout << "line-line intersection: " << p << endl;
    CGAL_assertion(p == sol);
    cout << "The solution is exact!" << endl;


    NT root24 = sqrt(sqrt(F2));
    NT root234 = root24*root24*root24;
    NT x1  = (NT(-15815)*sqrt(F2) + NT(19119)*root24 + NT(13444)*root234 - NT(23479))*sqrt(NT(1)+sqrt(F2))/NT(542) + NT(24831)*sqrt(F2)/NT(542) - NT(14840)*root24/NT(271) - NT(10708)*root234/NT(271) + NT(17856)/NT(271);
    NT y1  = (NT(-9090)*sqrt(F2) + NT(9803)*root24 + NT(8313)*root234 - NT(13801))*sqrt(NT(1)+sqrt(F2))/NT(542) + NT(14453)*sqrt(F2)/NT(542) - NT(8526)*root24/NT(271) - NT(6195)*root234/NT(271) + NT(9501)/NT(271);
    NT r12 = (NT(-3028175552)*sqrt(F2) + NT(3601721772)*root24 + NT(2545971914)*root234 - NT(4279338600))*sqrt(NT(1)+sqrt(F2))/NT(146882) + NT(4705259439)*sqrt(F2)/NT(146882) - NT(2797188994)*root24/NT(73441) - NT(1977303520)*root234/NT(73441) + NT(3327746973)/NT(73441);

    Circle c1(Point(x1, y1), r12);

    std::pair<Point, Point> ipt = Construct_intersection_2()(ell1, c1);
    Point ip1 = ipt.first, ip2 = ipt.second;
    cout << "Intersection of circle and line: " << ip1 << " and " << ip2 << endl;
    CGAL_assertion(ip1 == sol || ip2 == sol);
    cout << "The solution is exact!" << endl;


    NT x2  = (NT(203109)*sqrt(F2) - NT(54251)*root24 + NT(10397)*root234 + NT(223071))*sqrt(sqrt(F2) - NT(1))/NT(84386) - NT(48258)*sqrt(F2)/NT(42193) + NT(23571)*root24/NT(42193) + NT(106401)*root234/NT(84386) - NT(14511)/NT(42193);
    NT y2  = (NT(-406218)*sqrt(F2) + NT(277274)*root24 + NT(105785)*root234 - NT(446142))*sqrt(F2 - sqrt(F2))/NT(168772) - NT(23571)*sqrt(F2)/NT(42193) + NT(14511)*root24/NT(42193) + NT(48258)*root234/NT(42193) - NT(106401)/NT(42193);
    NT r22 = (NT(-38737365474)*sqrt(F2) + NT(62674657740)*root24 + NT(43541114565)*root234 - NT(69201833202))*sqrt(NT(1)+sqrt(F2))/NT(3560498498) + NT(87751517789)*sqrt(F2)/NT(3560498498) - NT(35806387098)*root24/NT(1780249249) - NT(28113144981)*root234/NT(1780249249) + NT(125172318975)/NT(3560498498);

    Circle c2(Point(x2, y2), r22);
    std::pair<Point, Point> cpt = Construct_intersection_2()(c1, c2);
    ip1 = cpt.first;
    ip2 = cpt.second;
    cout << "Intersection of circles: " << ip1 << " and " << ip2 << endl;
    CGAL_assertion(ip1 == sol || ip2 == sol);
    cout << "The solution is exact!" << endl;

    return 0;
}
