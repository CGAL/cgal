#include <iostream>
#include <boost/cstdint.hpp>
#include <boost/static_assert.hpp>

int main()
{
  std::cout << "Verifying the sizes of boost::[u]int{8,16,32,64}_t"
            << std::endl;
  
  BOOST_STATIC_ASSERT(sizeof(boost::int8_t)  == 1);
  BOOST_STATIC_ASSERT(sizeof(boost::int16_t) == 2);
  BOOST_STATIC_ASSERT(sizeof(boost::int32_t) == 4);
#ifndef BOOST_NO_INT64_T
  BOOST_STATIC_ASSERT(sizeof(boost::int64_t) == 8);
#endif

  BOOST_STATIC_ASSERT(sizeof(boost::uint8_t)  == 1);
  BOOST_STATIC_ASSERT(sizeof(boost::uint16_t) == 2);
  BOOST_STATIC_ASSERT(sizeof(boost::uint32_t) == 4);
#ifndef BOOST_NO_INT64_T
  BOOST_STATIC_ASSERT(sizeof(boost::uint64_t) == 8);
#endif

  return 0;
}
