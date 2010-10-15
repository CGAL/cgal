// Description: This program computes the miniball of some random
// balls in 3D using the LP-algorithm and double arithmetic, and thus
// checks that the LP-algorithm at least compiles and runs using
// double arithmetic (although this test doesn't say anything about
// the quality of the produced result).
//
// Note: This program is identical to the benchmark.cpp program in
// the examples directory, with the only difference that the LP-algorithm
// is run instead of the default algorithm.  (In particular, the output
// radii should be the same.)
//
// Compatibility: works with or without CGAL

#include <CGAL/basic.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

#include <iostream>
#include <vector>
#include <ctime>

class Ball {
private: // representation:
  double c[3]; // center in Eucliden coordinates
  double r;    // radius

public: // constructor:
  Ball() {}

  template<typename InputIterator>
  Ball(InputIterator from,double r) : r(r) {
    c[0] = *from;
    c[1] = *++from;
    c[2] = *++from;
  }

public: // accessors:
  double radius() const { return r; }

public: // iterator to iterate over the 3 coordinates:
  typedef const double *ConstIterator;
  ConstIterator beginCenter() const { return c; }
};

struct BallTraits {
  typedef ::Ball Sphere;
  static const int D=3;
  typedef CGAL::LP_algorithm Algorithm;
  typedef CGAL::Tag_true Use_square_roots;
  typedef double FT;

  typedef Sphere::ConstIterator Cartesian_const_iterator;
  static Cartesian_const_iterator center_cartesian_begin(const Sphere& b) {
    return b.beginCenter();
  }
  static double radius(const Sphere& b) {
    return b.radius();
  }
};

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
  const int N = 1000000;
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
  mb.is_empty();
  time = clock() - time;

  // output running time:
  cout << "----------------------------------------------------" << endl
       << "Benchmark: "
       << static_cast<double>(time)/CLOCKS_PER_SEC
       << "s (in units of " << 1.0/CLOCKS_PER_SEC << "s)." << endl
       << "----------------------------------------------------" << endl
       << endl;

  // for the fun of it, output the radius:
  cout << "Done.  (Radius: " << mb.radius() << ')' << endl;
}
