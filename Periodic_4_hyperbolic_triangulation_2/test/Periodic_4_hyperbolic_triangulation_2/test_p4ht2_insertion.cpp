
#define PROFILING_MODE

#include <CGAL/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_random_points_in_disc_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

#include <CGAL/Timer.h>

#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

typedef CORE::Expr                                                              NT;

typedef CGAL::Cartesian<NT>                                                     Kernel;

//typedef CGAL::Circular_kernel_2<K, CGAL::Algebraic_kernel_for_circles_2_2<K> >  Kernel;


typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                    CGAL::Hyperbolic_octagon_translation>       Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef CGAL::Hyperbolic_octagon_translation_matrix<NT>                         Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_original_octagon                                        Side_of_original_octagon;
typedef Triangulation::Face_iterator                                            Face_iterator;

typedef CGAL::Cartesian<double>::Point_2                                        Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                          Creator;


long calls_apply_identity(0);
long calls_apply_non_identity(0);
long calls_append_identity(0);
long calls_append_non_identity(0);

long calls_predicate_identity(0);
long calls_predicate_non_identity(0);
double time_predicate_identity(0);
double time_predicate_non_identity(0);
double time_remove_dp(0);

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

    Side_of_original_octagon pred;

    cout << "---- for best results, make sure that you have compiled me in Release mode ----" << endl;
    
    double extime = 0.0;

    for (int exec = 1; exec <= iters; exec++) {
        std::vector<Point_double> v;
        std::vector<Point> pts;
        // We can consider points only in the circle circumscribing the fundamental domain 
        Hyperbolic_random_points_in_disc_2_double(v, 5*N, -1, 0.159);

        int cnt = 0;
        int idx = 0;
        do {    
            Point pt = Point(v[idx].x(), v[idx].y());
            if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
                pts.push_back(pt);
                cnt++;
            } 
            idx++;
        } while (cnt < N && idx < v.size());

        if (cnt < N) {
            cout << "I failed at the generation of random points! Exiting..." << endl;
            return -1;
        }

        cout << "iteration " << exec << ": inserting into triangulation (rational dummy points)... "; cout.flush();
        Triangulation tr;
        tr.insert_dummy_points(true);  

        CGAL::Timer tt;
        tt.start();
        tr.insert(pts.begin(), pts.end(), true);
        tt.stop();
        cout << "DONE! (# of vertices = " << tr.number_of_vertices() << ", time = " << tt.time() << " secs)" << endl;
        extime += tt.time();

        int bfc = 0;
        for (Face_iterator fit = tr.faces_begin(); fit != tr.faces_end(); fit++) {
            if (!(fit->translation(0).is_identity() && fit->translation(1).is_identity() && fit->translation(2).is_identity())) {
                bfc++;
            }
        }
        int Nf = tr.number_of_faces();
        double perc = (double)bfc/(double)Nf*100.0;
        cout << "Total number of faces      : " << Nf << endl;
        cout << "Faces crossing the boundary: " << bfc << endl;
        cout << "Percentage                 : " << perc << endl;
    
            cout << "Triangulation is valid: " << (tr.is_valid() ? "YES" : "NO") << endl;
    }

    double avgtime = extime / (double)iters;
    cout << "---------------------------------------" << endl;
    cout << "Average execution time over " << iters << " iterations: " << avgtime << " secs" << endl << endl;

    
    cout << "Calls to append resulting in     identity: " << calls_append_identity << endl;
    cout << "Calls to append resulting in non-identity: " << calls_append_non_identity << endl;
    cout << "Percentage                               : " << ((double)(calls_append_non_identity)/(double)(calls_append_non_identity+calls_append_identity)*100.0) << endl << endl;
    cout << "Calls to apply  with             identity: " << calls_apply_identity << endl;
    cout << "Calls to apply  with         non-identity: " << calls_apply_non_identity << endl;
    cout << "Percentage                               : " << ((double)(calls_apply_non_identity)/(double)(calls_apply_non_identity+calls_apply_identity)*100.0) << endl << endl;

    cout << "Predicate calls with     identity translations: " << calls_predicate_identity << endl;
    cout << "Predicate calls with non-identity translations: " << calls_predicate_non_identity << endl;
    cout << "Percentage                               : " << (double)calls_predicate_non_identity/(double)(calls_predicate_non_identity+calls_predicate_identity)*100.0 << endl << endl;

    cout << "Time in predicates with     identity translations: " << time_predicate_identity << endl;
    cout << "Time in predicates with non-identity translations: " << time_predicate_non_identity << endl;
    cout << "Percentage                                  : " << (double)time_predicate_non_identity/(double)(time_predicate_non_identity+time_predicate_identity)*100.0 << endl << endl;

    cout << "Time to remove dummy points                 : " << time_remove_dp << endl;

    return 0;
}
