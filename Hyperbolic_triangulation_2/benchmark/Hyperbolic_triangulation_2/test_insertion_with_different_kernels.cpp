
#include <CGAL/basic.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Timer.h>


typedef CORE::Expr                                                    NT1;
typedef CGAL::Cartesian<NT1>                                          Kernel1;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<Kernel1>     Traits1;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits1>            Triangulation1;
typedef Triangulation1::Point                                         P1;
typedef Triangulation1::Face_handle                                   FH1;
typedef Triangulation1::Voronoi_point                                 VP1;

typedef CGAL::Exact_circular_kernel_2                                 Kernel2;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel2>  Traits2;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits2>            Triangulation2;
typedef Triangulation2::Point                                         P2;
typedef Triangulation2::Face_handle                                   FH2;
typedef Triangulation2::Voronoi_point                                 VP2;


typedef CGAL::Cartesian<double>             Dummy_K;
typedef Dummy_K::Point_2                    CP;
typedef CGAL::Creator_uniform_2<double,CP>  Creator;

using std::cout;
using std::endl;

int main(int argc, char** argv) {    

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_points_to_insert] [optional: number_of_iterations]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);

    int iterations = 1;
    if (argc > 2) {
        iterations = atoi(argv[2]);
    }

    std::vector<double> extCK, extCORE;

    for (int iter = 0; iter < iterations; iter++) {
        std::cout << "--> generating random points..."; std::cout.flush();
        std::vector<CP> pts;
        CGAL::Random_points_in_disc_2<CP,Creator> g( 1.0 );
        CGAL::cpp11::copy_n( g, N, std::back_inserter(pts));
        std::cout << "  DONE!" << std::endl;

        std::cout << "--> copying random points...   "; std::cout.flush();
        std::vector<P1> pts1;
        std::vector<P2> pts2;
        for (int i = 0; i < pts.size(); i++) {
            pts1.push_back( P1(pts[i].x(), pts[i].y()) );
            pts2.push_back( P2(pts[i].x(), pts[i].y()) );
        }
        std::cout << "  DONE!" << std::endl;

        CGAL::Timer tmr;
        tmr.reset();

        Triangulation1 tr1;
        std::cout << "--> inserting into triangulation with CK...        "; std::cout.flush();
        tmr.start();
        tr1.insert( pts1.begin(), pts1.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << tr1.number_of_vertices() << std::endl;
        double t1 = tmr.time();

        tmr.reset();

        Triangulation2 tr2;
        std::cout << "--> inserting into triangulation with CORE::Expr..."; std::cout.flush();
        tmr.start();
        tr2.insert( pts2.begin(), pts2.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << tr2.number_of_vertices() << std::endl;
        double t2 = tmr.time();

        std::cout << std::endl;
        std::cout << "Time for CK:         " << t1 << std::endl;
        std::cout << "Time for CORE::Expr: " << t2 << std::endl;
        std::cout << std::endl << "===========================================================" << std::endl << std::endl;
        extCK.push_back(t1);
        extCORE.push_back(t2);
    }

    double avgCK = 0;
    double avgCORE = 0;
    for (int i = 0; i < iterations; i++) {
        avgCK += extCK[i];
        avgCORE += extCORE[i];
    }
    avgCK = avgCK/(double)iterations;
    avgCORE = avgCORE/(double)iterations;

    std::cout << "######################################################################" << std::endl << std::endl;
    std::cout << "Average time for CK:   " << avgCK   << std::endl;
    std::cout << "Average time for CORE: " << avgCORE << std::endl;
    std::cout << std::endl;

    return 0;
}
