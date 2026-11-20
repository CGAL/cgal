#ifndef SDG_MAIN_H
#define SDG_MAIN_H

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>

#include "test.h"

int main()
{
  CGAL::test_x(std::cin, "bizarre", false);
  CGAL::test_no_x(std::cin, "bizarre", false);

  CGAL::test_x(std::cin, "sitesx", false);
  CGAL::test_no_x(std::cin, "sitesx", false);

  CGAL::test_x(std::cin, "sitesxx", false);
  CGAL::test_no_x(std::cin, "sitesxx", false);

  CGAL::test_x(std::cin, "MartinHeldBugreport", false);
  CGAL::test_no_x(std::cin, "MartinHeldBugreport", false);

#include "inp_sdg.h"

  return 0;
}

#endif // SDG_MAIN_H
