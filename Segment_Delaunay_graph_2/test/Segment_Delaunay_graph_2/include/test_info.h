#ifndef CGAL_TEST_INFO_H
#define CGAL_TEST_INFO_H 1

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Random.h>
#include "IO/io_aux.h"

namespace CGAL {

template<class SDG>
bool test_info(SDG& sdg, const char* fname)
{
  CGAL::Random r(static_cast<int>(0));

  std::ifstream ifs(fname);
  assert( ifs );

  sdg.clear();
  typename SDG::Site_2  site =
    SDG::Site_2::construct_site_2(typename SDG::Point_2(CGAL::ORIGIN));

  // read the sites and insert them in the segment Delaunay graph
  int info_id = 1;
  std::cout << "Input:" << std::endl;
  while ( ifs >> site ) {
    Multi_info<int> info = info_id;
    info_id++;
    std::cout << "SITE TO BE INSERTED: "
              << site << " " << info << std::endl;
    sdg.insert(site, info);
  }
  std::cout << std::endl;

  typedef typename SDG::Finite_vertices_iterator FVIT;
  for (FVIT it = sdg.finite_vertices_begin();
       it != sdg.finite_vertices_end(); ++it) {
    std::cout << it->site() << " "
              << it->storage_site().info() << std::endl;
  }
  std::cout << std::endl;

  // validate the segment Delaunay graph
  bool valid = sdg.is_valid(true, 1);
  std::cout << std::endl;
  return valid;
}

} //namespace CGAL

#endif // CGAL_TEST_INFO_H
