// Copyright (c) 2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Remy Thomasse  <remy.thomasse@inria.fr>

#ifndef CGAL_RANDOM_CONVEX_HULL_DISC_H
#define CGAL_RANDOM_CONVEX_HULL_DISC_H 1
#include <boost/random.hpp>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/function_objects.h>
#include <CGAL/copy_n.h>
#include <CGAL/number_type_config.h>
#include <list>

namespace CGAL {
namespace internal {
template <class P, class Traits>
struct compare_points_angle {
  bool operator()(const P& p, const P& q) {
    P zero(0,0);
    Traits traits;
    typedef typename Traits::Orientation_2 Orientation_2;
    Orientation_2 orientation_2=traits.orientation_2_object();
    typedef typename Traits::Compare_y_2 Compare_y_2;
    Compare_y_2 compare_y_2=traits.compare_y_2_object();

    if (compare_y_2(p, zero) == LARGER ){
      if (compare_y_2(q, zero)==LARGER)
        return (orientation_2(zero,p,q)==LEFT_TURN);
      else
        return false;

    } else {
      if (compare_y_2(q,zero)==LARGER)
        return true;
      else
        return (orientation_2(zero,p,q)==LEFT_TURN);
    }
  }
};

//////////////////////////////////////
template <class P, class GEN>
void generate_points_annulus(long n, double a, double b, double small_radius,
                             double big_radius, std::list<P>& l,
                             GEN& gen) {  // generate n points between a and b
  if (n > 1) {
    boost::binomial_distribution<long> bin_distribution(n, .5);
    boost::variate_generator<GEN&, boost::binomial_distribution<long> >
        g(gen, bin_distribution);
    long nb = g();
    generate_points_annulus(nb, a, (a + b) / 2.0, small_radius, big_radius, l,
                            gen);
    generate_points_annulus(n - nb, (a + b) / 2.0, b, small_radius, big_radius,
                            l, gen);
  }
  if (n == 1)  // generation of a point
  {

    #if BOOST_VERSION < 104700

    boost::uniform_real<double> random_squared_radius_distribution(
        small_radius * small_radius / (big_radius * big_radius), 1);
    boost::uniform_real<double> random_angle_distribution(a, b);
    boost::variate_generator<
        GEN&, boost::uniform_real<double> > random_angle(gen, random_angle_distribution);
    boost::variate_generator<
        GEN&, boost::uniform_real<double> > random_squared_radius(gen, random_squared_radius_distribution);

    #else

    boost::random::uniform_real_distribution<double> random_squared_radius_distribution(
        small_radius * small_radius / (big_radius * big_radius), 1);

    boost::random::uniform_real_distribution<double> random_angle_distribution(a, b);
    boost::random::variate_generator<
        GEN&, boost::random::uniform_real_distribution<double> > random_angle(gen, random_angle_distribution);
    boost::random::variate_generator<
        GEN&, boost::random::uniform_real_distribution<double> > random_squared_radius(gen, random_squared_radius_distribution);

    #endif

    double alpha = random_angle();
    double r = big_radius * std::sqrt(random_squared_radius());
    typedef  Creator_uniform_2<double, P> Creator;
    Creator creator;
    typedef typename Creator::argument_type T;
    l.push_back(creator(T(r * cos(alpha)), T(r * std::sin(alpha))));
  }
}

template <class P>
void Cyclic_increment(typename std::list<P>::iterator& it,
                               std::list<P>& l) {
  ++it;
  if (it == l.end()) {
    it = l.begin();
  }
}
//////////////////////////////////////////////////////////////////////////////
template <class P, class Traits>
void Graham_without_sort_2(std::list<P>& l, const Traits& traits) {
  if (l.size() > 3) {
     //typedef typename Traits::Left_turn_2 Left_turn;
     //Left_turn left_turn = traits.left_turn_2_object();
    typedef typename Traits::Orientation_2 Orientation_2;
    typedef typename std::list<P>::iterator Iterator;
    Orientation_2 orientation_2=traits.orientation_2_object();

     typedef typename Traits::Compare_x_2 Compare_x_2;
     Compare_x_2 compare_x_2=traits.compare_x_2_object();

    Iterator pmin = l.begin();
    for (Iterator it = l.begin(); it != l.end(); ++it) {
      if (compare_x_2(*pmin, *it) == LARGER){
        pmin = it;
      }
    }  //*pmin is the extremal point on the left
    Iterator u = pmin;
    Iterator u_next = u;
    Cyclic_increment(u_next, l);

    Iterator u_next_next = u_next;
    Cyclic_increment(u_next_next, l);

    while (u_next != pmin) {

      if (orientation_2(*u,*u_next,*u_next_next)==LEFT_TURN){
        Cyclic_increment(u, l);
        Cyclic_increment(u_next, l);
        Cyclic_increment(u_next_next, l);
      } else {
        u_next = l.erase(u_next);
        if (u_next == l.end()) u_next = l.begin();
        if (u != pmin) {
          u_next = u;
          if (u == l.begin()) {
            u = l.end();
          }
          --u;
        } else {
          u_next_next = u_next;
          Cyclic_increment(u_next_next, l);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
template <class GEN, class Traits>

void random_convex_hull_in_disc_2(std::size_t n, double radius, std::list<typename Traits::Point_2>& l,
                                  GEN& gen, const Traits& traits,
                                  bool fast = true) {
  CGAL_precondition(n >= 3);
  typedef typename Traits::Point_2 P;
  typedef typename Traits::FT FT;
  std::size_t simulated_points = 0;
  std::size_t generated_points = 0;
  P zero(0, 0);
  FT squared_radius = radius * radius;
  typedef typename Traits::Compare_y_2 Compare_y_2;
  Compare_y_2 compare_y_2=traits.compare_y_2_object();



  do {  // Initialization
    long init =
        static_cast<long>((std::min)(static_cast<std::size_t>(100), n - simulated_points));

    generate_points_annulus(init, -CGAL_PI, CGAL_PI, 0, radius, l,
                            gen);

    simulated_points += init;
    generated_points += init;

    Graham_without_sort_2(l, traits);
  } while ((bounded_side_2(l.begin(), l.end(), zero, traits) !=
            ON_BOUNDED_SIDE) &&
           (simulated_points < n));  // initialization such that 0 in P_n

  std::size_t T = n;
  if (!fast) T = static_cast<std::size_t>(std::floor(n / std::pow(std::log(static_cast<double>(n)), 2)));

  while (simulated_points < n) {
    // l is a list coming from a convex hull operation. we are moving the
    // points s.t the angles are from -pi to pi.
    {
      typename std::list<P>::iterator it = l.begin();
      while (compare_y_2(*it,zero) == LARGER){
        l.push_back(*it);
        l.pop_front();
        it = l.begin();
      }
      it = l.end();
      --it;  // last element
      while (compare_y_2(*it,zero) == SMALLER){
        l.push_front(*it);
        l.pop_back();
        it = l.end();
        --it;  // last element
      }
    }
    FT squared_small_radius = squared_radius;

    {
      typename std::list<P>::iterator it = l.begin();
      typename std::list<P>::iterator it2=++l.begin();
      for (; it != l.end();
           ++it, Cyclic_increment(it2, l)) {  // computation of annulus
        typename Traits::Segment_2 s(*it, *it2);
        FT temp=squared_distance(s,zero);
        if ( compare(squared_small_radius,temp) == LARGER ) squared_small_radius=temp;
      }
    }  // squared_small_radius=squared small radius of the annulus

    FT p_disc = squared_small_radius / squared_radius;
    long nb;
    if (simulated_points < T) {
      nb = static_cast<long>((std::min)(simulated_points, n - simulated_points));
    } else {
      nb = static_cast<long>((std::min)(T, n - simulated_points));
    }
    boost::binomial_distribution<long> dbin(nb, to_double(p_disc));
    boost::variate_generator<
        GEN&, boost::binomial_distribution<long> > bin(gen, dbin);

    // How many points are falling in the small disc and wont be generated:
    long k_disc = bin();
    simulated_points += k_disc;

    std::list<P> m;
    generate_points_annulus(nb - k_disc, -CGAL_PI, CGAL_PI,
                                      std::sqrt(to_double(squared_small_radius)), radius,
                                      m, gen);
    l.merge(m, compare_points_angle<P, Traits>());
    generated_points += nb - k_disc;
    simulated_points += nb - k_disc;
    m.clear();
    Graham_without_sort_2(l, traits);
  }
}

}  // namespace CGAL::internal

///

template <class OutputIterator, class Traits, class Generator>
void random_convex_hull_in_disc_2(std::size_t n, double radius, Generator& gen,
                                  OutputIterator it, const Traits& traits,
                                  bool fast = true) {
  typedef typename Traits::Point_2 Points;
  std::list<Points> l;
  internal::random_convex_hull_in_disc_2(n, radius, l, gen, traits, fast);
  std::copy_n(l.begin(),l.size(),it);
}

}  // namespace CGAL
#endif
