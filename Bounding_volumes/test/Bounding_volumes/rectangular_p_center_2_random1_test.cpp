// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA
#include <vector>
#include <functional>
#include <algorithm>
#include <cstdlib>
#ifndef CGAL_PCENTER_NO_OUTPUT
#include <iostream>
#include <iterator>
#endif // CGAL_PCENTER_NO_OUTPUT

using std::vector;
using std::back_inserter;
using std::atoi;
using std::exit;
using CGAL::Cartesian;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::Random;
using CGAL::Timer;
using CGAL::rectangular_p_center_2;
#ifndef CGAL_PCENTER_NO_OUTPUT
using std::ostream;
using std::cerr;
using std::flush;
using std::endl;
using std::copy;
using std::ostream_iterator;
#endif // CGAL_PCENTER_NO_OUTPUT

#ifdef CGAL_USE_LEDA
typedef leda_real                          FT;
#else
typedef double                             FT;
#endif // CGAL_USE_LEDA

typedef Cartesian< FT >                    K;
typedef K::Point_2                         Point;
typedef K::Vector_2                        Vector;
typedef K::Iso_rectangle_2                 Square_2;
typedef vector< Point >                    PCont;
typedef PCont::iterator                    iterator;
typedef Creator_uniform_2< FT, Point >     Creator;
typedef Random_points_in_square_2< Point, Creator >
                                           Point_generator;
template < class P,
           class Creator =
           CGAL::Creator_uniform_2< typename P::R::FT, P > >
class Random_p_clusters_2 : public CGAL::Random_generator_base< P > {
  void generate_point() {
    typedef typename P::R::FT FT;
    double p = this->_rnd.get_double();
    Creator creator;
    if (p <= 1.0 / n)
      d_item =
        creator(FT(p0.x() + c_size * (2 * this->_rnd.get_double() - 1.0)),
                FT(p0.y() + c_size * (2 * this->_rnd.get_double() - 1.0)));
    else if (p <= 2.0 / n)
      d_item =
        creator(FT(p1.x() + c_size * (2 * this->_rnd.get_double() - 1.0)),
                FT(p1.y() + c_size * (2 * this->_rnd.get_double() - 1.0)));
    else if (p <= 3.0 / n)
      d_item =
        creator(FT(p2.x() + c_size * (2 * this->_rnd.get_double() - 1.0)),
                FT(p2.y() + c_size * (2 * this->_rnd.get_double() - 1.0)));
    else
      d_item =
        creator(FT(p3.x() + c_size * (2 * this->_rnd.get_double() - 1.0)),
                FT(p3.y() + c_size * (2 * this->_rnd.get_double() - 1.0)));
  }

public:

  typedef Random_p_clusters_2< P, Creator > This;
  typedef CGAL::Random_generator_base< P >  Base;

  using Base::d_item;
  using Base::d_range;

  Random_p_clusters_2(int n_,
                      double c_size_,
                      double r = 1,
                      Random& rnd = CGAL::get_default_random())
  : Base(r - c_size_, rnd),
    n(n_),
    c_size(c_size_),
    p0(Creator()(d_range * (2 * this->_rnd.get_double() - 1.0),
                 d_range * (2 * this->_rnd.get_double() - 1.0))),
    p1(Creator()(d_range * (2 * this->_rnd.get_double() - 1.0),
                 d_range * (2 * this->_rnd.get_double() - 1.0))),
    p2(Creator()(d_range * (2 * this->_rnd.get_double() - 1.0),
                 d_range * (2 * this->_rnd.get_double() - 1.0))),
    p3(Creator()(d_range * (2 * this->_rnd.get_double() - 1.0),
                 d_range * (2 * this->_rnd.get_double() - 1.0)))
  {
    assert(n >= 1 && n <= 4);
    assert(c_size >= 0 && c_size <= r);
    generate_point();
  }
  This& operator++() {
    generate_point();
    return *this;
  }
  This  operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
private:
  int n;
  double c_size;
  P p0, p1, p2, p3;
};

int
main(int argc, char* argv[])
{
#ifndef CGAL_PCENTER_NO_OUTPUT
  CGAL::IO::set_pretty_mode(cerr);
#endif // CGAL_PCENTER_NO_OUTPUT

  int number_of_points;
  int random_seed;

  // handle command line arguments:
  if (argc < 2 || (number_of_points = atoi(argv[1])) <= 0) {
    cerr << "usage: " << argv[0]
         << " num [random_seed]" << endl;
    exit(1);
  }
  if (argc < 3) {

#ifndef CGAL_PCENTER_NO_OUTPUT
    cerr << "No random seed specified - generating it" << endl;
#endif // CGAL_PCENTER_NO_OUTPUT

    // generate random seed
    random_seed = CGAL::get_default_random().get_int(0, (1 << 30));
  }
  else
    random_seed = atoi(argv[2]);

  // define random source:
  Random rnd(random_seed);

#ifndef CGAL_PCENTER_NO_OUTPUT
  cerr << "random seed is " << random_seed << endl;
#endif // CGAL_PCENTER_NO_OUTPUT
  PCont input_points;
  std::copy_n(Point_generator(1, rnd),
                number_of_points,
                back_inserter(input_points));

  for (int p(2); p < 5; ++p) {
#ifndef CGAL_PCENTER_NO_OUTPUT
    cerr << "** computing " << p << "-centers:" << endl;
#endif // CGAL_PCENTER_NO_OUTPUT

    PCont centers;
    FT result;
    Timer t;
    t.start();
    rectangular_p_center_2(
      input_points.begin(),
      input_points.end(),
      back_inserter(centers),
      result,
      p);
    t.stop();
#ifndef CGAL_PCENTER_NO_OUTPUT
    cerr << "[time: " << t.time() << " sec]" << endl;
#endif // CGAL_PCENTER_NO_OUTPUT

#ifdef CGAL_USE_LEDA
    // check that all points are covered
    CGAL::I_Infinity_distance_2< K > dist;
    #ifdef _MSC_VER
    {
    #endif
    for (iterator i = input_points.begin(); i != input_points.end(); ++i) {
      iterator j = centers.begin();
      do {
        if (dist(*i, *j) <= result / FT(2))
          break;
        if (++j == centers.end()) {
    #ifndef _MSC_VER
          cerr << "Error: Point " << *i << " is not covered." << endl;
    #else
          cerr << "Error: A point is not covered." << endl;
    #endif
          assert(j != centers.end());
        }
      } while (j != centers.end());
    }
    #ifdef _MSC_VER
    }
    #endif

    // check that there is at least one square with two points
    // on opposite sides
    CGAL::I_Signed_x_distance_2< K > xdist;
    CGAL::I_Signed_y_distance_2< K > ydist;
    bool boundary = false;
    #ifdef _MSC_VER
    {
    #endif
    for (iterator i = centers.begin(); i != centers.end(); ++i) {
      int left = 0, right = 0, bottom = 0, top = 0;
      for (iterator j = input_points.begin(); j != input_points.end(); ++j) {
        if (xdist(*i, *j) == result / FT(2))
          ++left;
        if (xdist(*j, *i) == result / FT(2))
          ++right;
        if (ydist(*j, *i) == result / FT(2))
          ++top;
        if (ydist(*i, *j) == result / FT(2))
          ++bottom;
      }
      if ( (left > 0 && right > 0) || (top > 0 && bottom > 0) )
        boundary = true;
    }
    #ifdef _MSC_VER
    }
    #endif
    if (!boundary)
      cerr << "Error: No square has two points on boundary." << endl;
    assert(boundary);
#endif // CGAL_USE_LEDA

  } // for (int p(2); p < 5; ++p)

  return 0;
}
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

