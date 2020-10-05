#include <boost/config.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <iostream>
#include <fstream>

#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <vector>
#include <deque>

//typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;


int main(int argc, char* argv[])
{
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef std::deque<Polygon> MultiPolygon;

  {
    std::ifstream is((argc>1)?argv[1]:"data/polygons.wkt");
    std::list<Polygon> polys;
    do
      {
        Polygon p;
        CGAL::read_polygon_WKT(is, p);
        if(!p.outer_boundary().is_empty())
          polys.push_back(p);
      }while(is.good() && !is.eof());
    for(Polygon p : polys)
      std::cout<<p<<std::endl;
  }

  {
    std::ifstream  is((argc>2)?argv[2]:"data/multipolygon.wkt");
    MultiPolygon mp;
    CGAL::read_multi_polygon_WKT(is, mp);
    for(Polygon p : mp)
      std::cout<<p<<std::endl;
  }
  return 0;
}
#else
int main()
{
  return 0;
}
#endif
