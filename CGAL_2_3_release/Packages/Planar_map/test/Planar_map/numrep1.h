// Including headers for testing
#if TESTR == 1 || TESTR == 3 || TESTR == 5 || TESTR == 7 
#include <CGAL/Cartesian.h>
#endif
#if TESTR == 2 || TESTR == 4 || TESTR == 6 || TESTR == 8
#include <CGAL/Homogeneous.h>
#endif

#if TESTR == 2
#if !defined(CGAL_USE_LEDA)
#define VOID_TEST
#define VOID_TEST_VENDOR "LEDA"
#else
#include <CGAL/leda_integer.h>
//#include <CGAL/leda_debug_integer.h>
#endif
#endif

#if TESTR == 3 || TESTR == 4
#if !defined(CGAL_USE_GMP)
#define VOID_TEST
#define VOID_TEST_VENDOR "GMPZ"
#else
#include <CGAL/Gmpz.h>
#endif
#endif

#if TESTR == 5 || TESTR == 6
#if !defined(CGAL_USE_LEDA)
#define VOID_TEST
#define VOID_TEST_VENDOR "LEDA"
#else
#include <CGAL/leda_rational.h>
#endif
#endif

#if TESTR == 7 || TESTR == 8
#include <CGAL/test_types.h>
#endif

#if TESTR == 2 || TESTR == 5 || TESTR == 6

#endif

#if TESTR == 1
inline double to_nt(int d)
{
    return (double) d;
}
#endif

#if TESTR == 2
inline double to_nt(int d)
{
    return d;
}
#endif

#if TESTR == 7
inline CGAL::TestfieldC to_nt(int d)
{
    unsigned char dummy1 = 'a';
    signed char dummy2 = 'a';
    return CGAL::TestfieldC(dummy1, dummy2, (double)d);
}
#endif

#if TESTR == 8
inline CGAL::TestrepH to_nt(int d)
{
    unsigned char dummy1 = 'a';
    signed char dummy2 = 'a';
    return CGAL::TestrepH(dummy1, dummy2, d);
}
#endif

struct randomint {
    randomint() ;
    int	get() const { return sequence[cur]; }
    int next() { cur = (cur+1)%11; return get();}
private:
    int sequence[11];
    int cur;
};

inline randomint::randomint()
{
    cur = 0;
    sequence[0] = 19;
    sequence[1] = 5;
    sequence[2] = 17;
    sequence[3] = 13;
    sequence[4] = 29;
    sequence[5] = 2;
    sequence[6] = 23;
    sequence[7] = 31;
    sequence[8] = 3;
    sequence[9] = 37;
    sequence[10] = 11;
}

#ifdef VOID_TEST
void usage(){
  std::cout << VOID_TEST_VENDOR << " installation missing.\n";
  std::cout << "Please install " << VOID_TEST_VENDOR << " for this configuration to be tested.\n";
  std::cout << "Test is not performed.\n";
}
#endif







