#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Object.h>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>               Traits;
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_with_history_2< Traits, Dcel> Arrangement;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement>
  Walk_along_line_pl;
typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Ray_2                                   Ray_2;
typedef Arrangement::Curve_2                            Curve_2;

int main( )
{
    Point_2 p1(0, 0);
    Point_2 p2(1, 0);
    Ray_2 ray(p1, p2);
    Curve_2 curve(ray);

    Arrangement arr;
    CGAL::insert(arr, curve);

    Walk_along_line_pl pl(arr);
    auto o = pl.locate(Point_2(1, -1));

    return 0;
}
