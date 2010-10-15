#include <iostream>
#include <ctime>
#include <CGAL/Min_sphere_of_spheres_d.h>

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
  typedef const double *Coord_iterator;
  Coord_iterator begin_center() const { return c; }
};

struct Ball_traits {
  typedef Ball Sphere;
  static const int D=3;
  typedef double FT;
  typedef CGAL::Default_algorithm Algorithm;
  typedef CGAL::Tag_false Use_square_roots;
  typedef Sphere::Coord_iterator Cartesian_const_iterator;

  static Cartesian_const_iterator
    center_cartesian_begin(const Ball& b) {
    return b.begin_center();
  }
  static double radius(const Ball& b) {
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
  typedef CGAL::Min_sphere_of_spheres_d<Ball_traits> Minsphere;
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
       << "Benchmark: " << static_cast<double>(time)/CLOCKS_PER_SEC
       << "s (in units of " << 1.0/CLOCKS_PER_SEC << "s)." << endl
       << "----------------------------------------------------" << endl
       << endl;

  // output the radius:
  cout << "Done.  (Radius: " << mb.radius() << ')' << endl;
}
