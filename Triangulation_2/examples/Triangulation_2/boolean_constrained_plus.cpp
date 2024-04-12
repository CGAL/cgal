#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>


#include <CGAL/Triangulation_2/Boolean.h>
#include <CGAL/draw_constrained_triangulation_2.h>
#include <CGAL/draw_multipolygon_with_holes_2.h>

#include <CGAL/IO/WKT.h>


#include <iostream>
#include <sstream>

#include <boost/property_map/property_map.hpp>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<K>;

int
main( )
{
  CGAL::Triangulations::Boolean<K> bops;

  Multipolygon_with_holes_2 pA, pB;

  {
    std::istringstream is("MULTIPOLYGON( ((0 0,  20 0, 20 30, 0 30, 0 0), (1 1, 1 2, 2 2, 2 1, 1 1 ) ) ,   (( 50 0, 60 0, 60 60, 50 60)) )");
    //std::istringstream is("MULTIPOLYGON( ((0 0, 2 0, 2 3, 0 3) ) )");    // (0.1 0.1, 0.1 0.4, 0.4 0.1)
    CGAL::IO::read_multi_polygon_WKT(is, pA);
  }

  {
    std::istringstream is("MULTIPOLYGON( ((10 1,  30 1, 30 2, 20 2, 20 4, 10 4)) )");
    //std::istringstream is("MULTIPOLYGON( ((2 1, 3 1, 3 2, 2 2)) ");
    CGAL::IO::read_multi_polygon_WKT(is, pB);
  }


  bops.insert(pA,pB);
  Multipolygon_with_holes_2 mpwh = bops([](bool a, bool b){ return a || b;});
  CGAL::IO::write_multi_polygon_WKT(std::cout, mpwh);
  CGAL::draw(mpwh);

  /*
    std::map<Boolean_cdt_2::Face_handle,bool> map;
    for(auto fh : bops.cdt.finite_face_handles()){
      map[fh] = fh->info().in_domain(0) || fh->info().in_domain(1);
      assert(map[fh] == (fh->info().label != 0));
   }

    CGAL::draw(bops.cdt, boost::make_assoc_property_map(map));
  */
  return 0;
}
