#include <cassert>
#include <iostream>

int main()
{
  std::cout << "Verifying the sizes of boost::[u]int{8,16,32,64}_t"
            << std::endl;

  assert(sizeof(boost::int8_t)  == 1);
  assert(sizeof(boost::int16_t) == 2);
  assert(sizeof(boost::int32_t) == 4);
#ifndef BOOST_NO_INT64_T
  assert(sizeof(boost::int64_t) == 8);
#endif

  assert(sizeof(boost::uint8_t)  == 1);
  assert(sizeof(boost::uint16_t) == 2);
  assert(sizeof(boost::uint32_t) == 4);
#ifndef BOOST_NO_INT64_T
  assert(sizeof(boost::uint64_t) == 8);
#endif

  return 0;
}
