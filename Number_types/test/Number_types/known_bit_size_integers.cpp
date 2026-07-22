#include <iostream>
#include <cstdint>
#include <CGAL/assertions.h>

int main()
{
  std::cout << "Verifying the sizes of std::[u]int{8,16,32,64}_t"
            << std::endl;

  static_assert(sizeof(std::int8_t)  == 1);
  static_assert(sizeof(std::int16_t) == 2);
  static_assert(sizeof(std::int32_t) == 4);
  static_assert(sizeof(std::int64_t) == 8);


  static_assert(sizeof(std::uint8_t)  == 1);
  static_assert(sizeof(std::uint16_t) == 2);
  static_assert(sizeof(std::uint32_t) == 4);
  static_assert(sizeof(std::uint64_t) == 8);

  return 0;
}
