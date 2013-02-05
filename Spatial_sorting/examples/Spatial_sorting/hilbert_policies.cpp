#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hilbert_sort.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;

int main ()
{
  std::vector<Point> v; v.reserve(4);
  v.push_back( Point(0.0,0.0)) ;
  v.push_back( Point(1.0,1.0)) ;
  v.push_back( Point(0.1,0.1)) ;
  v.push_back( Point(0.2,0.8)) ;
  
  std::cout << "Hilbert sort (middle policy)." << std::endl;
  CGAL::hilbert_sort (v.begin(), v.end(), K(), CGAL::Hilbert_sort_middle_policy());
  std::cout<<v[0]<<"; "<<v[1]<<"; "<<v[2]<<"; "<<v[3]<<"; "<<std::endl;
  std::cout << "Hilbert sort (median policy)." << std::endl;
  CGAL::hilbert_sort (v.begin(), v.end(), K(), CGAL::Hilbert_sort_median_policy());
  std::cout<<v[0]<<"; "<<v[1]<<"; "<<v[2]<<"; "<<v[3]<<"; "<<std::endl;
  return 0;
}
