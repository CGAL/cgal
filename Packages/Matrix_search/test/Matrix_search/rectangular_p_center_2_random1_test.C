// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : rectangular_p_center_2_random1_test.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Centering Axis-Parallel 2D-Rectangles - test program
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_VECTOR_2_H
#include <CGAL/Vector_2.h>
#endif // CGAL_VECTOR_2_H
#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RECTANGULAR_P_CENTER_2_H
#include <CGAL/rectangular_p_center_2.h>
#endif // CGAL_RECTANGULAR_P_CENTER_2_H
#ifndef CGAL_RANDOM_H
#include <CGAL/Random.h>
#endif // CGAL_RANDOM_H
#ifndef CGAL_ALGORITHM_H
#include <CGAL/algorithm.h>
#endif // CGAL_ALGORITHM_H
#include <vector>
#include <functional>
#include <algorithm>
#include <cstdlib>
#ifndef CGAL_PCENTER_NO_OUTPUT
#include <iostream>
#include <iterator>
#endif // CGAL_PCENTER_NO_OUTPUT

#include <ctime>
static time_t Measure;
static long long int measure;
#define MEASURE(comm) \
Measure = clock(); \
  comm; \
  measure = (long long int)((float)(clock() - Measure) \
  * 1000 / CLOCKS_PER_SEC); \
    cerr << "[time: " << measure << " msec]\n";
#define MEASURE_NO_OUTPUT(comm) \
Measure = clock(); \
  comm; \
  measure = (long long int)((float)(clock() - Measure) \
  * 1000 / CLOCKS_PER_SEC);

using std::vector;
//    using std::bind2nd;
//    using std::transform;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::Random;
using CGAL::default_random;
using CGAL::rectangular_p_center_2;
#ifndef CGAL_PCENTER_NO_OUTPUT
using std::ostream;
using std::cerr;
using std::flush;
using std::endl;
using std::copy;
using std::ostream_iterator;
#endif // CGAL_PCENTER_NO_OUTPUT

typedef double                             FT;
//#include <CGAL/Arithmetic_filter.h>
//#include <CGAL/leda_real.h>
//using CGAL::Filtered_exact;
//typedef Filtered_exact< double, leda_real >  FT;

typedef Cartesian< FT >                    R;
typedef CGAL::Point_2< R >                 Point;
typedef CGAL::Vector_2< R >                Vector;
typedef CGAL::Iso_rectangle_2< R >         Square_2;
typedef vector< Point >                    PCont;
typedef PCont::iterator                    iterator;
typedef Creator_uniform_2< FT, Point >     Creator;
typedef Random_points_in_square_2< Point, Creator >
                                           Point_generator;
#ifdef _MSC_VER
// that compiler cannot even distinguish between global
// and class scope, so ...
#define Base B_B_Base
#endif // _MSC_VER

#ifdef CGAL_CFG_NO_NAMESPACE
template < class P,
           class Creator =
           CGAL_STD::Creator_uniform_2< typename P::FT, P > >
class Random_p_clusters_2 : public CGAL_STD::_Random_generator_base< P > {
#else
template < class P,
           class Creator =
           CGAL::Creator_uniform_2< typename P::FT, P > >
class Random_p_clusters_2 : public CGAL::_Random_generator_base< P > {
#endif
  void generate_point() {
    typedef typename P::FT FT;
    double p = _rnd.get_double();
    Creator creator;
    if (p <= 1.0 / n)
      d_item =
        creator(FT(p0.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p0.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else if (p <= 2.0 / n)
      d_item =
        creator(FT(p1.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p1.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else if (p <= 3.0 / n)
      d_item =
        creator(FT(p2.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p2.y() + c_size * (2 * _rnd.get_double() - 1.0)));
    else
      d_item =
        creator(FT(p3.x() + c_size * (2 * _rnd.get_double() - 1.0)),
                FT(p3.y() + c_size * (2 * _rnd.get_double() - 1.0)));
  }
public:
  typedef Random_p_clusters_2< P, Creator > This;
  typedef CGAL::_Random_generator_base< P > Base;
  Random_p_clusters_2(int n_,
                      double c_size_,
                      double r = 1,
                      Random& rnd = default_random)
  : Base(r - c_size_, rnd),
    n(n_),
    c_size(c_size_),
    p0(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p1(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p2(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0))),
    p3(Creator()(d_range * (2 * _rnd.get_double() - 1.0),
                 d_range * (2 * _rnd.get_double() - 1.0)))
  {
    CGAL_precondition(n >= 1 && n <= 4);
    CGAL_precondition(c_size >= 0 && c_size <= r);
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

#ifdef _MSC_VER
#undef Base
#endif // _MSC_VER

int
main(int argc, char* argv[])
{
#ifndef CGAL_PCENTER_NO_OUTPUT
  CGAL::set_pretty_mode(cerr);
  ostream_iterator< Point > cerr_it_p(cerr, "\n");
  //ostream_iterator< Square_2 > cerr_it_s(cerr, "\n");
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
  cerr << "warning: no random seed specified\n"
         << "generating random seed" << endl;
#endif // CGAL_PCENTER_NO_OUTPUT

    // generate random seed
    random_seed = default_random.get_int(0, (1 << 30));
  }
  else
    random_seed = atoi(argv[2]);

  // define random source:
  Random rnd(random_seed);

#ifndef CGAL_PCENTER_NO_OUTPUT
  cerr << "random seed is " << random_seed << endl;
#endif // CGAL_PCENTER_NO_OUTPUT
  PCont input_points;
  CGAL::copy_n(Point_generator(1, rnd),
                number_of_points,
                back_inserter(input_points));

  for (int p(2); p <= 4; ++p) {
#ifndef CGAL_PCENTER_NO_OUTPUT
    cerr << "** computing " << p << "-centers:" << endl;
#endif // CGAL_PCENTER_NO_OUTPUT

    PCont centers;
    FT result;
    MEASURE(
      rectangular_p_center_2(
        input_points.begin(),
        input_points.end(),
        back_inserter(centers),
        result,
        p);
      )

    CGAL::Infinity_distance_2< R > dist;
    for (iterator i = input_points.begin(); i != input_points.end(); ++i) {
      iterator j = centers.begin();
      do {
        if (dist(*i, *j) <= result / FT(2))
          break;
        if (++j == centers.end()) {
          cerr << "Error: Point " << *i << " is not covered." << endl;
          CGAL_assertion(j != centers.end());
        }
      } while (j != centers.end());
    }

  } // for (int p(2); p < 4; ++p)

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

