#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/random_convex_hull_in_disc_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <boost/random.hpp>
#include <iostream>
#include <vector>
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel          K;
typedef K::Point_2                                       Point;
typedef K::FT                                                                                 FT;

const double RADIUS=1.0;
int main( )
{
   int N=10000;
   std::vector<Point> v;
   boost::mt19937 gen;
   gen.seed(0u);

   random_convex_hull_in_disc_2(N,RADIUS,gen,std::back_inserter(v),K());
   size_t size = v.size();
   FT area=polygon_area_2(v.begin(),v.end(),K());
   std::cout<<"A random convex polygon inscribed in a disc with "<<size<<" vertices and area "<<area<<" has been generated."<<std::endl;

   return 0;
}
