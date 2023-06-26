#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Random.h>


#include <sstream>

template <class K>
void test()
{
  typedef typename K::FT FT;
  CGAL::Delaunay_triangulation_3<K> dt3;


  // create a bunch of points
  CGAL::Random random;
  FT x,y,z;
  for (int n=0;n<50;++n)
  {
    x=random.get_int(-500,500);
    y=random.get_int(-500,500);
    z=random.get_int(-500,500);
    dt3.insert(typename K::Point_3(x/FT(3),y/FT(5),z/FT(7)));
  }

  std::stringstream buffer;
  buffer << std::setprecision(17) << dt3;
  decltype(dt3) dt3_bis;
  buffer >> dt3_bis;

  assert(dt3==dt3_bis);
}

int main()
{
  test<CGAL::Simple_cartesian<CGAL::Exact_rational>>();
  test<CGAL::Exact_predicates_exact_constructions_kernel>();
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();
}