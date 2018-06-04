
#include <CGAL/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Timer.h>

typedef double                                                                  NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>                                  Triangulation;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef CGAL::Creator_uniform_2<double, Point>                                  Creator;


typedef CORE::Expr                                                              NT2;
typedef CGAL::Cartesian<NT2>                                                    Kernel2;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel2,
                                      CGAL::Hyperbolic_octagon_translation>     Traits2;
typedef Traits2::Side_of_original_octagon                                       Side_of_original_octagon;

using std::cout;
using std::endl;

int main(int argc, char** argv) {    

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_points_to_insert] [optional: number_of_iterations]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);
    int iters = 1;
    if (argc == 3) {
        iters = atoi(argv[2]);
    }

    cout << "---- for best results, make sure that you have compiled me in Release mode ----" << endl;
    
    double extime = 0.0;
    Side_of_original_octagon pred;

    for (int exec = 1; exec <= iters; exec++) {
        std::vector<Point> pts;
        CGAL::Random_points_in_disc_2<Point, Creator> g(0.85);

        int cnt = 0;
        do {    
            Point pt = *(++g);
            if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
                pts.push_back(pt);
                cnt++;
            } 
        } while (cnt < N);

        cout << "iteration " << exec << ": inserting into triangulation (rational dummy points)... "; cout.flush();
        Triangulation tr;

        CGAL::Timer tt;
        tt.start();
        tr.insert(pts.begin(), pts.end());
        tt.stop();
        cout << "DONE! (# of vertices = " << tr.number_of_vertices() << ", time = " << tt.time() << " secs)" << endl;
        extime += tt.time();
    }

    double avgtime = extime / (double)iters;
    cout << "---------------------------------------" << endl;
    cout << "Average execution time over " << iters << " iterations: " << avgtime << " secs" << endl << endl;

    return 0;
}
