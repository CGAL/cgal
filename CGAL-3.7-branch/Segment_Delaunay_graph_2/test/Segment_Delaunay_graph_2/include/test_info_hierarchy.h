#ifndef CGAL_TEST_INFO_HIERARCHY_H
#define CGAL_TEST_INFO_HIERARCHY_H 1

#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include "IO/io_aux.h"
#include "test_info.h"

namespace CGAL {

template<class SDG>
bool test_info_hierarchy(SDG& sdg, const char* fname)
{ 
  typedef typename SDG::Finite_vertices_iterator   FVIT;

  bool valid = test_info(sdg, fname);

  std::cout << std::endl;
  std::cout << "ALL SITES/ALL LEVELS:" << std::endl;
  for (int i = 0; i <= 4; ++i) {
    std::cout << "\tSITES FOR LEVEL " << i << std::endl;
    for (FVIT it = sdg.diagram(i).finite_vertices_begin();
	 it != sdg.diagram(i).finite_vertices_end(); ++it) {
      std::cout << "\t\t" << it->site() << " "
		<< it->storage_site().info() << std::endl;
    }
  }
  std::cout << std::endl;

  return valid;
}

} //namespace CGAL

#endif // CGAL_TEST_INFO_HIERARCHY_H
