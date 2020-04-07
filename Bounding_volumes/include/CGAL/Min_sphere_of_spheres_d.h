// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Kaspar Fischer


#ifndef CGAL_MINIBALL_MINIBALL
#define CGAL_MINIBALL_MINIBALL

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_configure.h>
#include <boost/random/linear_congruential.hpp>
#include <cmath>
#include <vector>
#include <iostream>

#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_pair.h>
#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_support_set.h>

namespace CGAL_MINIBALL_NAMESPACE {

  // We provide several algorithms to solve the Miniball problem.  The
  // following types are used by the client to select the algorithm he
  // or she wants to run.
  struct LP_algorithm {};
  struct Farthest_first_heuristic {};
  typedef Farthest_first_heuristic Default_algorithm;

  template<class Traits>
  class Min_sphere_of_spheres_d {
  public: // some short hands:
    typedef typename Traits::Sphere Sphere;
    typedef typename Traits::FT FT;
    typedef typename Min_sphere_of_spheres_d_impl::
                     Selector<FT>::Result Result;
    typedef typename Min_sphere_of_spheres_d_impl::
                     Selector<FT>::Is_exact Is_exact;
    typedef typename Traits::Use_square_roots Use_sqrt;
    typedef typename Traits::Algorithm Algorithm;
    static const int D = Traits::D;
    typedef typename Traits::Cartesian_const_iterator CIt;

  private: // traits class:
    Traits t; // To allow the traits to not only vary at compile- but
              // also at runtime, we instantiate it here.

  private: // for internal consisteny checks:
    // The following variable is true if and only if the miniball
    // has been computed of all inserted balls, i.e. iff every checked-in
    // ball has been respected in the miniball computation.
    bool is_up_to_date;

    boost::rand48 rng;

  public: // iterators:
    typedef const Result *Cartesian_const_iterator; // coordinate iterator

    class Support_iterator {
      typedef typename std::vector<const Sphere*>::const_iterator It;
      It it;

    private:
      friend class Min_sphere_of_spheres_d<Traits>;
      Support_iterator(It it) : it(it) {}

    public:
      const Sphere& operator*() { return *(*it); }
      Support_iterator& operator++() { ++it; return *this; }
      Support_iterator operator++(int) {
        Support_iterator old(*this);
        ++(*this);
        return old;
      }
      bool operator!=(const Support_iterator& i) { return it != i.it; }
    };

  public: // construction and destruction:
    inline Min_sphere_of_spheres_d(const Traits& traits = Traits());

    template<typename InputIterator>
    inline Min_sphere_of_spheres_d(InputIterator begin,InputIterator end,
                                   const Traits& traits = Traits()) :
      t(traits), ss(t), e(0) {
      is_up_to_date = false;
      CGAL_MINIBALL_ASSERT(is_neg(ss.radius(),ss.disc()));
      insert(begin,end);            // todo. better way?
    }


    inline void insert(const Sphere& b);

    template<typename InputIterator>
    inline void insert(InputIterator begin,InputIterator end) {
      S.insert(S.end(), begin, end);
      is_up_to_date = false;
#ifdef CGAL_MINIBALL_DEBUG
      for(std::size_t i = 0; i < S.size(); i++){
        CGAL_MINIBALL_ASSERT(t.radius(S[i]) >= FT(0));
      }
#endif
    }

    inline void clear();

    template<typename InputIterator>
    inline void set(InputIterator begin,InputIterator end) {
      clear();
      insert(begin,end);
    }

  public: // predicates and accessors:
    inline bool is_empty();

    inline const Result& radius();

    inline Cartesian_const_iterator center_cartesian_begin();
    inline Cartesian_const_iterator center_cartesian_end();

    inline const FT& discriminant();

    inline Support_iterator support_begin();

    inline Support_iterator support_end();

    inline const Traits& traits() const {
      return t;
    }

  public: // validity check:
    bool is_valid();

  private:
    bool is_valid(const Tag_true is_exact);
    bool is_valid(const Tag_false is_exact);

  private:
    void update();
    void update(LP_algorithm);
    void update(Farthest_first_heuristic);

    bool find_farthest(int from,int to,int& i,
                       Tag_true use_sqrt,Tag_false is_exact);
    bool find_farthest(int from,int to,int& i,
                       Tag_true use_sqrt,Tag_true is_exact);
    bool find_farthest(int from,int to,int& i,
                       Tag_false use_sqrt,Tag_false is_exact);
    bool find_farthest(int from,int to,int& i,
                       Tag_false use_sqrt,Tag_true is_exact);

  private:
    std::vector<Sphere> S;         // list of the added bals
    std::vector<const Sphere *> l; // list of pointers to the added bals
    Min_sphere_of_spheres_d_impl::
      Support_set<Traits> ss;      // current support set
    int e;                         // l[0..(e-1)] is a basis

  private: // forbid copy constructing and assignment (because our
           // pointers in ss would be wrong then):
    Min_sphere_of_spheres_d(const Min_sphere_of_spheres_d&);
    Min_sphere_of_spheres_d& operator=(const Min_sphere_of_spheres_d&);
  };

  template<class Traits>
  Min_sphere_of_spheres_d<Traits>::
    Min_sphere_of_spheres_d(const Traits& traits) :
    t(traits), ss(t), e(0) {
    is_up_to_date = true;
    CGAL_MINIBALL_ASSERT(is_neg(ss.radius(),ss.disc())); // makes sure
               // that initially no ball is enclosed (cf. contains()).
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::insert(const Sphere& b) {
    CGAL_MINIBALL_ASSERT(t.radius(b) >= FT(0));
    S.push_back(b);
    // (We push_back a pointer to S.back() in update().)
    is_up_to_date = false;
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::clear() {
    S.clear();
    l.clear();
    ss.reset();
    e = 0;
    is_up_to_date = true;
  }

  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::is_empty() {
    if (!is_up_to_date)
      update();
    return is_neg(ss.radius(),ss.disc());
  }

  template<class Traits>
  const typename Min_sphere_of_spheres_d<Traits>::Result&
    Min_sphere_of_spheres_d<Traits>::radius() {
    if (!is_up_to_date)
      update();
    CGAL_MINIBALL_ASSERT(!is_empty());
    return ss.radius();
  }

  template<class Traits>
  typename Min_sphere_of_spheres_d<Traits>::Cartesian_const_iterator
    Min_sphere_of_spheres_d<Traits>::center_cartesian_begin() {
    if (!is_up_to_date)
      update();
    CGAL_MINIBALL_ASSERT(!is_empty());
    return ss.begin();
  }

  template<class Traits>
  typename Min_sphere_of_spheres_d<Traits>::Cartesian_const_iterator
    Min_sphere_of_spheres_d<Traits>::center_cartesian_end() {
    if (!is_up_to_date)
      update();
    CGAL_MINIBALL_ASSERT(!is_empty());
    return ss.begin()+D;
  }

  template<class Traits>
  const typename Min_sphere_of_spheres_d<Traits>::FT&
    Min_sphere_of_spheres_d<Traits>::discriminant() {
    if (!is_up_to_date)
      update();
    CGAL_MINIBALL_ASSERT(!is_empty());
    return ss.disc();
  }

  template<class Traits>
  inline void Min_sphere_of_spheres_d<Traits>::update() {
    // set up the vector l containing pointers to the balls in S:
    CGAL_MINIBALL_ASSERT(l.size() == 0);
    for (typename std::vector<Sphere>::const_iterator it = S.begin();
         it != S.end(); ++it)
      l.push_back(&(*it));

    // compute the miniball:
    update(Algorithm());
    is_up_to_date = true;
  }

  template<class Traits>
  inline bool Min_sphere_of_spheres_d<Traits>::is_valid() {
    if (!is_up_to_date)
      update();
    return is_valid(Is_exact());
  }

  template<class Traits>
  inline typename Min_sphere_of_spheres_d<Traits>::Support_iterator
    Min_sphere_of_spheres_d<Traits>::support_begin() {
    if (!is_up_to_date)
      update();
    return Support_iterator(l.begin());
  }

  template<class Traits>
  inline typename Min_sphere_of_spheres_d<Traits>::Support_iterator
    Min_sphere_of_spheres_d<Traits>::support_end() {
    if (!is_up_to_date)
      update();
    return Support_iterator(l.begin()+e);
  }

} // namespace CGAL_MINIBALL_NAMESPACE

// If the package is used with CGAL, we include some default
// traits classes:
#ifndef CGAL_JUST_MINIBALL
#include <CGAL/Min_sphere_of_spheres_d_traits_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#endif

#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_impl.h>
#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_pivot_impl.h>

#endif // CGAL_MINIBALL_MINIBALL
