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
using Vector_2 = K::Vector_2;
using Transformation = K::Aff_transformation_2;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<K>;

int
main(int argc, char* argv[])
{
  CGAL::Triangulations::Boolean<K> bops;

  Multipolygon_with_holes_2 pA, pB;
  if(argc == 2) {
      std::ifstream in(argv[1]);
      CGAL::IO::read_multi_polygon_WKT(in, pA);

      CGAL::Bbox_2 bb = pA.bbox();
      double w = bb.xmax() - bb.xmin();
      Vector_2 vec(w / 10.0, w / 10.0);
      pB = CGAL::transform(Transformation(CGAL::TRANSLATION, vec), pA);
  }
  else if (argc == 3) {
      {
          std::ifstream in(argv[1]);
          CGAL::IO::read_multi_polygon_WKT(in, pA);
      }
      {
          std::ifstream in(argv[2]);
          CGAL::IO::read_multi_polygon_WKT(in, pB);
      }
  }
  else {
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
  }

  bops.insert(pA,pB);
  Multipolygon_with_holes_2 mpwh = bops([](bool a, bool b){ return a || b;});
  std::ofstream out("result.wkt");
  CGAL::IO::write_multi_polygon_WKT(out, mpwh);
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
