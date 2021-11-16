#ifndef CUBE_H
#define CUBE_H

#include <array>
#include <cassert>
#include <string>
#include <sstream>

typedef std::array<unsigned char, 8> Cube;

inline constexpr Cube convert_to_cube(unsigned int n) {
  assert(n < (1<<24));
  return {
    (unsigned char)((n & 070000000)>>21), (unsigned char)((n & 007000000)>>18),
    (unsigned char)((n & 000700000)>>15), (unsigned char)((n & 000070000)>>12),
    (unsigned char)((n & 000007000)>> 9), (unsigned char)((n & 000000700)>> 6),
    (unsigned char)((n & 000000070)>> 3), (unsigned char)((n & 000000007)>> 0),
  };
}

// User-defined literal operator.  Given an integer in octal notation, like
// 01234567, gives the cube with the same colors. For example, `01234567_c`
// is `Cube{0, 1, 2, 3, 4, 5, 6, 7}`.
inline constexpr Cube operator"" _c ( unsigned long long n )
{
  assert(n < (1<<24));
  return convert_to_cube(unsigned(n));
}

inline std::string config_name(const Cube cube) {
  std::stringstream filename_prefix_ss;
  for(int j = 0; j < 8; ++j) {
    filename_prefix_ss << int(cube[j]);
  }
  return filename_prefix_ss.str();
}

#endif // CUBE_H
