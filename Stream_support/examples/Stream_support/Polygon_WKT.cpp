#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Multipolygon_with_holes_2.h>

#include <iostream>
#include <fstream>
#include <deque>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

int main(int argc, char* argv[])
{
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef std::deque<Polygon> MultiPolygon;
  typedef CGAL::Multipolygon_with_holes_2<Kernel> Multipolygon_with_holes_2;
  {
    std::ifstream is((argc>1)?argv[1]:"data/polygons.wkt");
    std::list<Polygon> polys;
    do
      {
        Polygon p;
        CGAL::IO::read_polygon_WKT(is, p);
        if(!p.outer_boundary().is_empty())
          polys.push_back(p);
      }while(is.good() && !is.eof());
    for(Polygon p : polys)
      std::cout<<p<<std::endl;
  }

  {
    std::ifstream  is((argc>2)?argv[2]:"data/multipolygon.wkt");
    MultiPolygon mp;
    CGAL::IO::read_multi_polygon_WKT(is, mp);
    for(Polygon p : mp)
      std::cout<<p<<std::endl;
  }

  {
    std::ifstream  is((argc>2)?argv[2]:"data/multipolygon.wkt");
    Multipolygon_with_holes_2 mp;
    CGAL::IO::read_multi_polygon_WKT(is, mp);
    std::cout << mp << std::endl;
    CGAL::IO::write_multi_polygon_WKT(std::cout, mp);
    std::cout << std::endl;
  }
  return 0;
}
