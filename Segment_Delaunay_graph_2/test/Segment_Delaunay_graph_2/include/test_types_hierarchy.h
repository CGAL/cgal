#ifndef CGAL_SDG_TEST_TYPES_HIERARCHY_H
#define CGAL_SDG_TEST_TYPES_HIERARCHY_H

#include "test_types.h"

namespace CGAL {

template<class SDG, class InputStream>
bool test_sdg_hierarchy(InputStream& is, const SDG& sdg, char* fname)
{
  return test_sdg(is, sdg, fname);
}

} //namespace CGAL

#endif // CGAL_SDG_TEST_TYPES_HIERARCHY_H
