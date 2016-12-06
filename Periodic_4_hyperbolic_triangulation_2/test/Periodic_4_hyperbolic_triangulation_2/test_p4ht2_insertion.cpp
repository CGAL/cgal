

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_random_points_in_disc_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

#include <time.h>



typedef CORE::Expr                                                              NT;         
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Hyperbolic_octagon_translation_matrix<Traits>                           Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_fundamental_octagon                                     Side_of_fundamental_octagon;

typedef CGAL::Cartesian<double>::Point_2                                        Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                          Creator;


int main(void) {    

    Triangulation tr_rational;//, tr_algebraic;    
    
    Side_of_fundamental_octagon pred;

    int N = 10000;

    tr_rational.insert_dummy_points(true);  

    std::vector<Point> pts;
    
    vector<Point_double> v;
    Hyperbolic_random_points_in_disc_2_double(v, 40*N, -1);

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
        return -1;
    }


    clock_t t_start_r = clock();
    for (int j = 0; j < N; j++) {
        Vertex_handle vh = tr_rational.insert(pts[j]);
    }
            
    
    double s_t1 = (double)(clock() - t_start_r)/CLOCKS_PER_SEC;
    

    cout << "Insertion time for " << N << " points:  " << s_t1 << endl;
    cout << "Vertices in triangulation: " << tr_rational.number_of_vertices() << endl;    

    return 0;
}


