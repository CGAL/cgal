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
// file          : include/CGAL/Min_sphere_of_spheres_d.h
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


#ifndef CGAL_MINIBALL_MINIBALL
#define CGAL_MINIBALL_MINIBALL

#include <CGAL/Min_sphere_of_spheres_d_configure.h>
#include <cassert>

#include <vector>
#include <iostream>
#include <CGAL/Min_sphere_of_spheres_d_support_set.h>

namespace CGAL_MINIBALL_NAMESPACE {

  struct LP_algorithm {};
  struct Farthest_first_heuristic {};
  typedef Farthest_first_heuristic Default_algorithm;

  template<typename FT>
  inline bool compare(const FT& a,const FT& b,
                      const FT& ap,const FT& bp) {
    const FT u = a-ap, uu = u*u;
    if (u >= FT(0)) {
      if (bp <= uu)
        return false;
  
      // here (1) holds
      const FT v = uu-b+bp;
      if (v <= 0)
        return false;
  
      // here (2) holds
      return 4 * uu * bp < sqr(v);
    } else {
      // here (1) holds
      const FT v = uu-b+bp;
      if (v >= FT(0))
        return true;
  
      // here (3) holds
      return 4 * uu *bp > sqr(v);
    }
  }

  template<class Traits>
  class Min_sphere_of_spheres_d {
  public: // some short hands:
    typedef typename Traits::Sphere Sphere;
    typedef typename Traits::FT FT;
    typedef typename Selector<FT>::Result Result;
    typedef typename Selector<FT>::Is_exact Is_exact;
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
  
  public: // iterators:
    typedef const Result *Cartesian_const_iterator; // coordinate iterator
  
    class Support_iterator {
      typedef typename std::vector<Sphere>::const_iterator It;
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
      t(traits), e(0), ss(t) {
      is_up_to_date = false;
      CGAL_MINIBALL_ASSERT(is_neg(ss.radius(),ss.disc()));
      insert(begin,end);            // todo. better way?
    }
  
    inline void prepare(int size);
  
    inline void insert(const Sphere& b);
  
    template<typename InputIterator>
    inline void insert(InputIterator begin,InputIterator end) {
      prepare(l.size()+(end-begin)); // todo. istream?
      while (begin != end) {
        insert(*begin);
        ++begin;
      }
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
    bool is_valid(const Tag_true is_exact);
    bool is_valid(const Tag_false is_exact);
  
  private:
    bool pivot(int B);
  
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
    std::vector<Sphere> l;        // list of the added bals
    Support_set<Traits> ss;       // current support set
    int e;                        // l[0..(e-1)] is a basis
  
  private: // forbid copy constructing and assignment (because our
           // pointers in ss would be wrong then):
    Min_sphere_of_spheres_d(const Min_sphere_of_spheres_d&);
    Min_sphere_of_spheres_d& operator=(const Min_sphere_of_spheres_d&);
  };

  template<class Traits>
  Min_sphere_of_spheres_d<Traits>::
    Min_sphere_of_spheres_d(const Traits& traits) :
    t(traits), e(0), ss(t) {
    is_up_to_date = true;
    CGAL_MINIBALL_ASSERT(is_neg(ss.radius(),ss.disc())); // makes sure
               // that initially no ball is enclosed (cf. contains()).
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::prepare(int size) {
    l.reserve(size);
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::insert(const Sphere& b) {
    CGAL_MINIBALL_ASSERT(t.radius(b) >= FT(0));
    l.push_back(b);
    is_up_to_date = false;
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::clear() {
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
#ifdef CGAL_VERSION
#include <CGAL/Min_sphere_of_spheres_d_traits_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#endif

#ifdef CGAL_MINIBALL_NO_TEMPLATE_EXPORT
#include <CGAL/Min_sphere_of_spheres_d.C>
#include <CGAL/Min_sphere_of_spheres_d_pivot.C>
#endif

#endif // CGAL_MINIBALL_MINIBALL
