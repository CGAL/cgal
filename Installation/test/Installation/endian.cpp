// Run-time check that our endianness macro is correctly defined.

#include <CGAL/config.h>
#include <cassert>

#if !defined CGAL_LITTLE_ENDIAN && !defined CGAL_BIG_ENDIAN
#  error no endian macro defined
#endif

union T {
  int       testWord;
  char      testByte[sizeof(int)];
} endianTest;

int main()
{
    endianTest.testWord = 1;
    if (endianTest.testByte[0] == 1) {
#ifdef CGAL_LITTLE_ENDIAN
        return 0;
#else
        assert(false);
#endif
    } else {
#ifdef CGAL_BIG_ENDIAN
        return 0;
#else
        assert(false);
#endif
    }
    return 0;
}
