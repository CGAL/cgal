#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef K::Point_2                                          Point;
typedef CGAL::Creator_uniform_2<double,Point>               Creator_2;

void compute_delaunay(std::vector<Point>::iterator it,
                        std::vector<Point>::iterator e){
    DT dt;
    DT::Face_handle hint;
    for( ;it!=e; ++it)  hint = dt.insert(*it, hint)->face();
    return;
}

void test_orders(std::vector<Point>::iterator b,
                        std::vector<Point>::iterator e){
    CGAL::Timer cost;

    cost.reset();cost.start();
    std::cout << "  Delaunay without spatial sort...                "<< std::flush;
    compute_delaunay(b,e);cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;

    cost.reset();cost.start();
    std::cout << "  Delaunay with Hilbert sort (median policy)...   " << std::flush;
    CGAL::hilbert_sort(b,e, CGAL::Hilbert_sort_median_policy());
    compute_delaunay(b,e);cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;

    cost.reset();cost.start();
    std::cout << "  Delaunay with Hilbert sort (middle policy)...   " << std::flush;
    CGAL::hilbert_sort(b,e, CGAL::Hilbert_sort_middle_policy() );
    compute_delaunay(b,e);cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;

    cost.reset();cost.start();
    std::cout << "  Delaunay with spatial sort (median policy)...   " << std::flush;
    CGAL::spatial_sort(b,e, CGAL::Hilbert_sort_median_policy());
    compute_delaunay(b,e);cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;

    cost.reset();cost.start();
    std::cout << "  Delaunay with spatial sort (middle policy)...   " << std::flush;
    CGAL::spatial_sort(b,e, CGAL::Hilbert_sort_middle_policy() );
    compute_delaunay(b,e);cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;
}

int main ()
{   int size1 = 1000,size2=100000;
    std::vector<Point> v,w;
    v.reserve(size1);w.reserve(size2);

    std::cout <<size1<< " points on a parabola" << std::endl;
    for (int i=0; i< size1; ++i) {
      double x= -size1 +i;
      v.push_back( Point( x, x*x ));
    }
    test_orders(v.begin(), v.end());

    CGAL::Random random (42);
    CGAL::Random_points_in_square_2<Point> gen (1.0, random);
    std::cout <<size2<< " points at random" << std::endl;
    for (int i=0; i< size2; ++i) {
      w.push_back( *gen++ );
    }
    test_orders(w.begin(), w.end());
}
