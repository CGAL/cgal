// Description: This program computes the miniball of some random
// balls in 3D using the LP-algorithm and exact arithmetic.  At the
// end of the computation, the program checks that the computed ball
// is indeed the smallest enclosing ball.  So this program serves as a
// test for the implementation of the LP-algorithm under exact
// arithmetic.
//
// Compatibility: works with or without CGAL

#include <iostream>
#include <ctime>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include "Rational.h"

// The program will work with numbers of type FT.
// You can change FT to leda_real, for instance.
typedef Rational FT;

class Ball {
private: // representation:
  FT c[3]; // center in Eucliden coordinates
  FT r;    // radius

public: // constructor:
  Ball() {}

  template<typename InputIterator>
  Ball(InputIterator from,const FT& r) : r(r) {
    c[0] = *from;
    c[1] = *++from;
    c[2] = *++from;
  }

public: // accessors:
  const FT& radius() const { return r; }

public: // iterator to iterate over the 3 coordinates:
  typedef const FT *ConstIterator;
  ConstIterator beginCenter() const { return c; }
};

struct BallTraits {
  typedef Ball Sphere;
  static const int D=3;
  typedef CGAL::LP_algorithm Algorithm;
  typedef CGAL::Tag_false Use_square_roots;
  typedef ::FT FT;

  typedef Sphere::ConstIterator Coordinate_iterator;
  static Coordinate_iterator begin(const Sphere& b) {
    return b.beginCenter();
  }
  static const FT& radius(const Sphere& b) {
    return b.radius();
  }
};

double to_double(const std::pair<FT,FT>& p,const FT& disc) {
  return to_double(p.first)+to_double(p.second)*
    std::sqrt(to_double(disc));
}

double uniform() {  // a (platform independent) random number generator
  static int lastNo = 230575L;
  const int a = 16807L, m = 2147483647L, q = 127773L, r = 2836L;
  int gamma = a * (lastNo % q) - r * (lastNo / q);
  if (gamma > 0)
    lastNo = gamma;
  else
    lastNo = gamma + m;
  return static_cast<double>(lastNo)/m;
}

int main(int,char **) {
  typedef CGAL::Min_sphere_of_spheres_d<BallTraits> Minsphere;
  using namespace std;

  // generate a million random spheres:
  const int N = 1000;
  vector<Ball> S;
  for (int i=0; i<N; ++i) {
    double a[3] = { uniform(), uniform(), uniform() };
    S.push_back(Ball(a,uniform()));
  }

  // remember current time:
  clock_t time = clock();

  // check in the balls:
  Minsphere mb(S.begin(),S.end());
  
  // measure time:
  bool e = mb.is_empty();
  time = clock() - time;
  
  // output running time:
  cout << "----------------------------------------------------" << endl
       << "Benchmark: "
       << static_cast<double>(time)/CLOCKS_PER_SEC
       << "s (in units of " << 1.0/CLOCKS_PER_SEC << "s)." << endl
       << "----------------------------------------------------" << endl
       << endl;
  
  // for the fun of it, output the radius:
  Minsphere::Result radius = mb.radius();
  cout << "(Radius: " << to_double(radius,mb.discriminant()) << ')' << endl;

  // check whether the computed ball is indeed the minsphere:
  bool valid = mb.is_valid();
  if (valid) 
    cout << "Validity check passed." << endl;
  else 
    cout << "WARNING: Validity check *not* passed!" << endl
	 << "Please contact the author <fischerk@inf.ethz.ch>!" << endl;
  return valid? 0 : 1;
}
