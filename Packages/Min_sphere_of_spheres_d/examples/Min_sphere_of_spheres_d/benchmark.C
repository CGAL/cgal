// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind.
//
// Every use of CGAL requires a license.
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation.
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : examples/Min_sphere_of_spheres_d/benchmark.C
// package       : Min_sphere_of_spheres_d (1.10)
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Kaspar Fischer
// maintainer    : Kaspar Fischer <fischerk@inf.ethz.ch>
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


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
  typedef Sphere::Coord_iterator Coordinate_iterator;

  static Coordinate_iterator begin(const Ball& b) {
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

  // compute the miniball:
  mb.update();

  // measure time:
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
