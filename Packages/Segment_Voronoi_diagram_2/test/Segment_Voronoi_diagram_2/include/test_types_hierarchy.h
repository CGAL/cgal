#ifndef CGAL_SVD_TEST_TYPES_HIERARCHY_H
#define CGAL_SVD_TEST_TYPES_HIERARCHY_H

#include "test_types.h"

CGAL_BEGIN_NAMESPACE

template<class SVD, class InputStream>
bool test_svd_hierarchy(InputStream& is, const SVD& svd)
{
  return test_svd(is, svd);
}

CGAL_END_NAMESPACE

#endif // CGAL_SVD_TEST_TYPES_HIERARCHY_H
