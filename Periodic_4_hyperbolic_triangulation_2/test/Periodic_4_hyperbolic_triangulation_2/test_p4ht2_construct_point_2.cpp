
#include <CGAL/basic.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>


typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                      CGAL::Hyperbolic_octagon_translation>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Kernel::Point_2                                                         Point;

typedef Traits::Construct_hyperbolic_point_2                                    CP2;
typedef Traits::Hyperbolic_translation                                          Trans;

using namespace std;
using namespace CGAL;

int main(int argc, char** argv) {    

    std::vector<Trans> gens;
    Trans::generators(gens);

    Point O(0,0);
    for(int i=0; i<gens.size(); ++i) {
        Point pt = CP2()(O, gens[i]);
        cout << "Image " << i << ": " << pt << endl;
    }

    return 0;
}
