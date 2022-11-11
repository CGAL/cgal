#include <iostream>
#include <boost/cstdint.hpp>
#include <CGAL/assertions.h>

int main()
{
  std::cout << "Verifying the sizes of boost::[u]int{8,16,32,64}_t"
            << std::endl;

  CGAL_static_assertion(sizeof(boost::int8_t)  == 1);
  CGAL_static_assertion(sizeof(boost::int16_t) == 2);
  CGAL_static_assertion(sizeof(boost::int32_t) == 4);
#ifndef BOOST_NO_INT64_T
  CGAL_static_assertion(sizeof(boost::int64_t) == 8);
#endif

  CGAL_static_assertion(sizeof(boost::uint8_t)  == 1);
  CGAL_static_assertion(sizeof(boost::uint16_t) == 2);
  CGAL_static_assertion(sizeof(boost::uint32_t) == 4);
#ifndef BOOST_NO_INT64_T
  CGAL_static_assertion(sizeof(boost::uint64_t) == 8);
#endif

  return 0;
}
