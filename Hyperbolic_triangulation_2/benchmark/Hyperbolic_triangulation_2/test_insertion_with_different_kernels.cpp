
#include <CGAL/basic.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Timer.h>


typedef CORE::Expr                                                     CORE_NumberType;
typedef CGAL::Cartesian<CORE_NumberType>                               CORE_Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<CORE_Kernel>  CORE_Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<CORE_Traits>         CORE_Triangulation;
typedef CORE_Triangulation::Point                                      CORE_Point;

typedef CGAL::Exact_circular_kernel_2                                  CK_Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<CK_Kernel> CK_Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<CK_Traits>           CK_Triangulation;
typedef CK_Triangulation::Point                                        CK_Point;


typedef CGAL::Delaunay_triangulation_2<CORE_Kernel>                    CORE_Euclidean_Triangulation;
//typedef CORE_Euclidean_Triangulation::Point_2                          CORE_Euclidean_Point;

typedef CGAL::Delaunay_triangulation_2<CK_Kernel>                      CK_Euclidean_Triangulation;
//typedef CK_Euclidean_Triangulation::Point_2                            CK_Euclidean_Point;

typedef CGAL::Cartesian<double>                                        Dummy_Kernel;
typedef Dummy_Kernel::Point_2                                          Double_Point;
typedef CGAL::Creator_uniform_2<double,Double_Point>                   Creator;

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

    std::vector<double> CK_Time, CORE_Time, CK_Euclidean_Time, CORE_Euclidean_Time;

    for (int iter = 0; iter < iterations; iter++) {
        std::cout << "--> generating random points..."; std::cout.flush();
        std::vector<Double_Point> pts;
        CGAL::Random_points_in_disc_2<Double_Point,Creator> g( 1.0 );
        CGAL::cpp11::copy_n( g, N, std::back_inserter(pts));
        std::cout << "  DONE!" << std::endl;

        std::cout << "--> copying random points...   "; std::cout.flush();
        std::vector<CORE_Point>           Vector_CORE;
        std::vector<CK_Point>             Vector_CK;
        //std::vector<CK_Euclidean_Point>   Vector_Euclidean_CK;
        //std::vector<CORE_Euclidean_Point> Vector_Euclidean_CORE;
        for (int i = 0; i < pts.size(); i++) {
            Vector_CK.push_back( CK_Point(pts[i].x(), pts[i].y()) );
            Vector_CORE.push_back( CORE_Point(pts[i].x(), pts[i].y()) );
            //Vector_Euclidean_CK.push_back( CK_Euclidean_Point(pts[i].x(), pts[i].y()) );
            //Vector_Euclidean_CORE.push_back( CORE_Euclidean_Point(pts[i].x(), pts[i].y()) );
        }
        std::cout << "  DONE!" << std::endl;

        CGAL::Timer tmr;

        tmr.reset();
        CK_Euclidean_Triangulation euclidean_CK_triangulation;
        std::cout << "--> inserting into Euclidean triangulation with CK...         "; std::cout.flush();
        tmr.start();
        euclidean_CK_triangulation.insert( Vector_CK.begin(), Vector_CK.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << euclidean_CK_triangulation.number_of_vertices() << std::endl;
        double Time_Euclidean_CK = tmr.time();

        tmr.reset();
        CK_Triangulation CK_triangulation;
        std::cout << "--> inserting into hyperbolic triangulation with CK...        "; std::cout.flush();
        tmr.start();
        CK_triangulation.insert( Vector_CK.begin(), Vector_CK.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << CK_triangulation.number_of_vertices() << std::endl;
        double Time_CK = tmr.time();

        tmr.reset();
        CORE_Euclidean_Triangulation euclidean_CORE_triangulation;
        std::cout << "--> inserting into Euclidean triangulation with CORE::Expr... "; std::cout.flush();
        tmr.start();
        euclidean_CORE_triangulation.insert( Vector_CORE.begin(), Vector_CORE.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << euclidean_CORE_triangulation.number_of_vertices() << std::endl;
        double Time_Euclidean_CORE = tmr.time();

        tmr.reset();
        CORE_Triangulation CORE_triangulation;
        std::cout << "--> inserting into hyperbolic triangulation with CORE::Expr..."; std::cout.flush();
        tmr.start();
        CORE_triangulation.insert( Vector_CORE.begin(), Vector_CORE.end() );
        tmr.stop();
        std::cout << "  DONE! Triangulation vertices: " << CORE_triangulation.number_of_vertices() << std::endl;
        double Time_CORE = tmr.time();

        std::cout << std::endl;
        std::cout << "Time for Euclidean CK:    " << Time_Euclidean_CK   << std::endl;
        std::cout << "Time for hyperbolic CK:   " << Time_CK             << std::endl;
        std::cout << "Time for Euclidean CORE:  " << Time_Euclidean_CORE << std::endl;
        std::cout << "Time for hyperbolic CORE: " << Time_CORE           << std::endl;
        std::cout << std::endl << "===========================================================" << std::endl << std::endl;
        CK_Time.push_back(Time_CK);
        CORE_Time.push_back(Time_CORE);
        CK_Euclidean_Time.push_back(Time_Euclidean_CK);
        CORE_Euclidean_Time.push_back(Time_Euclidean_CORE);
    }

    double avgCK = 0;
    double avgCORE = 0;
    double avgECK = 0;
    double avgECORE = 0;
    for (int i = 0; i < iterations; i++) {
        avgCK += CK_Time[i];
        avgCORE += CORE_Time[i];
        avgECK += CK_Euclidean_Time[i];
        avgECORE += CORE_Euclidean_Time[i];
    }
    avgCK    = avgCK/(double)iterations;
    avgCORE  = avgCORE/(double)iterations;
    avgECK   = avgECK/(double)iterations;
    avgECORE = avgECORE/(double)iterations;

    std::cout << "######################################################################" << std::endl << std::endl;
    std::cout << "Average time for Euclidean with CK:    " << avgECK   << std::endl;
    std::cout << "Average time for hyperbolic with CK:   " << avgCK    << std::endl;
    std::cout << "Average time for Euclidean with CORE:  " << avgECORE << std::endl;
    std::cout << "Average time for hyperbolic with CORE: " << avgCORE  << std::endl;
    std::cout << std::endl;

    return 0;
}
