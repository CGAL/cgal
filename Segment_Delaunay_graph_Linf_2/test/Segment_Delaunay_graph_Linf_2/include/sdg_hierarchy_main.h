#ifndef SDG_MAIN_HIERARCHY_H
#define SDG_MAIN_HIERARCHY_H

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>

#include "test.h"

int main(int, char**)
{
  CGAL::test_hierarchy_x(std::cin, "bizarre", false);
  CGAL::test_hierarchy_no_x(std::cin, "bizarre", false);

  CGAL::test_hierarchy_x(std::cin, "sitesx", false);
  CGAL::test_hierarchy_no_x(std::cin, "sitesx", false);

  CGAL::test_hierarchy_x(std::cin, "sitesxx", false);
  CGAL::test_hierarchy_no_x(std::cin, "sitesxx", false);

  CGAL::test_hierarchy_x(std::cin, "MartinHeldBugreport", false);
  CGAL::test_hierarchy_no_x(std::cin, "MartinHeldBugreport", false);

#include "inp_hier.h"

  return 0;
}

#endif // SDG_MAIN_HIERARCHY_H
