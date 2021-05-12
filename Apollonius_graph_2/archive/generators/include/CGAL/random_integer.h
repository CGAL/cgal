#ifndef CGAL_RANDOM_INTEGER_H
#define CGAL_RANDOM_INTEGER_H

#include <CGAL/basic.h>
#include <cstdlib>
#include <CGAL/Gmpz.h>

// type "man {rand, random, drand48}" for C functions that produce
// random numbers

//extern "C" int getpid();
//extern "C" int srandom(unsigned);
//extern "C" long random();

namespace CGAL {

// powers of 2 from 2^0 to 2^53
double
P2[54]={1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0,
        512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0, 32768.0,
        65536.0, 131072.0, 262144.0, 524288.0, 1048576.0,
        2097152.0, 4194304.0, 8388608.0, 16777216.0, 33554432.0,
        67108864.0, 134217728.0, 268435456.0, 536870912.0,
        1073741824.0, 2147483648.0, 4294967296.0, 8589934592.0,
        17179869184.0, 34359738368.0, 68719476736.0,
        137438953472.0, 274877906944.0, 549755813888.0,
        1099511627776.0, 2199023255552.0, 4398046511104.0,
        8796093022208.0, 17592186044416.0, 35184372088832.0,
        70368744177664.0, 140737488355328.0, 281474976710656.0,
        562949953421312.0, 1125899906842624.0, 2251799813685248.0,
        4503599627370496.0, 9007199254740992.0};

// set the random number generator seed
void set_seed(unsigned int seed)
{
  srandom(seed);
}

// return random integer of b bits
double random_integer(int b, bool allow_negative = true)
{
  double value;

  if (b > 27) {
    value =  ((int)random()) % ((int)P2[27]);
    value += (  ((int)random()) % ((int)P2[b - 27])  ) * P2[27];
  } else {
    value = ((int)random()) % (int)P2[b];
  }
  if ( allow_negative ) {
    value *= ( ((int) random()) % 2 == 0 ) ? 1 : -1;
  }
  return value;
}

template<class Random>
double random_integer(Random& r, unsigned int b, bool allow_negative = true)
{
  // returns random integers in the range [0, 2^b - 1).
  // b is required to be at least 1 and at most 52
  // and if allow_negative is true then the range includes negative
  // numbers as well and becomes: [-2^b + 1, 2^b - 1).
  CGAL_precondition( b >= 0 && b <= 52 );

  if ( b == 0 ) { return 0; }

  double M = pow(2.0,b);
  CGAL::Gmpz z;

  if ( allow_negative ) {
    z = r.get_double(-M, M);
  } else {
    z = r.get_double(0, M);
  }
  return CGAL::to_double(z);
}


template<class Random>
double random_even_integer(Random& r, unsigned int b,
                           bool allow_negative = true)
{
  // returns random even integers in the range [0, 2^b - 1).
  // b is required to be at least 1 and at most 52
  // and if allow_negative is true then the range includes negative
  // numbers as well and becomes: [-2^b + 1, 2^b - 1).
  CGAL_precondition( b >= 0 && b <= 52 );

  if ( b == 0 ) { return 0; }

  double M = pow(2.0,b);
  CGAL::Gmpz z;
  CGAL::Gmpz two(2);

  do {
    if ( allow_negative ) {
      z = r.get_double(-M, M);
    } else {
      z = r.get_double(0, M);
    }
  } while ( z % two != 0 );

  return CGAL::to_double(z);
}


} //namespace CGAL


#endif // CGAL_RANDOM_INTEGER_H
