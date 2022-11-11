#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
void compute_delaunay(std::vector<K::Point_2>::iterator it,
                        std::vector<K::Point_2>::iterator e){
    DT dt;
    DT::Face_handle hint;
    for( ;it!=e; ++it)  hint = dt.insert(*it, hint)->face();
}
int main ()
{   int size = 1000;
    std::vector<K::Point_2> v;
    v.reserve(size);
    CGAL::Timer cost;
    std::cout <<size<< " points on a parabola" << std::endl;
    for (int i=0; i< size; ++i) {
      double x= -size +i;
      v.push_back( K::Point_2( x, x*x ));
    }
    cost.reset();cost.start();
    std::cout << "  Delaunay without spatial sort... "<< std::flush;
    compute_delaunay(v.begin(),v.end());cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;
    cost.reset();cost.start();
    std::cout << "  Delaunay with Hilbert sort...    " << std::flush;
    CGAL::hilbert_sort(v.begin(),v.end());
    compute_delaunay(v.begin(),v.end());cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;
    cost.reset();cost.start();
    std::cout << "  Delaunay with spatial sort...    " << std::flush;
    CGAL::spatial_sort(v.begin(),v.end());
    compute_delaunay(v.begin(),v.end());cost.stop();
    std::cout << "done in "<<cost.time()<<" seconds." << std::endl;
    return 0;
}
