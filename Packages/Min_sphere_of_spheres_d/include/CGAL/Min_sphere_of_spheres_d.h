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
// release       : CGAL-2.4
// release_date  : 2002, May 16
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres.h
// package       : Min_sphere_of_spheres_d (1.00)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Kaspar Fischer
//
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_D_H
#define CGAL_MIN_SPHERE_OF_SPHERES_D_H

#include <cassert>
#include <vector>
#include <list>

#include <CGAL/basic.h>
#include <CGAL/Min_sphere_of_spheres_support_set.h>

namespace CGAL {

  // The following types select which algorithm is run, i.e. whether
  // the LP-algorithm or the "farthest-first" heuristic is used to
  // find the miniball.

  struct LP_algorithm {};
  struct Farthest_first_heuristic {};
  typedef Farthest_first_heuristic Default_algorithm;
  
  // another namespace to "hide" implementation details:
  namespace Min_sphere_of_spheres_impl {

    // The following routines are used to normalize rational numbers
    // during exact computation.  The assume that FT comes with a
    // Number_type_traits traits-class, check whether gcd() is
    // available, and if so, normalize the number:

    template<typename FT>
    void normalize(FT& a) {
      // a.normalize();
    }

    template<typename FT>
    void normalize(Pair<FT>& p) {
      normalize(p.first);
      normalize(p.second);
    }

    void normalize(double) {
    }

    void normalize(Pair<double>&) {
    }
    
    template<typename FT>
    inline bool compare(const FT& a,const FT& b,
			const FT& ap,const FT& bp) {
      const FT u = a-ap, uu = u*u;
      if (u >= 0) {
	if (bp <= uu)
	  return false;
	
	// here (1) holds
	const FT v = b-bp;
	if (uu - v <= 0)
	  return false;
	
	// here (2) holds
	return 4 * uu * bp < CGAL_NTS square(uu-v);
      } else {
	// here (1) holds
	const FT v = b-bp;
	if (uu-v >= 0)
	  return true;
	
	// here (3) holds
	return 4 * uu *bp > CGAL_NTS square(uu-v);
      }
    }

  } // namespace Min_sphere_of_spheres_impl
    
  template<class Traits>
  class Min_sphere_of_spheres_d {
  public: // public types:
    typedef typename Traits::Sphere Sphere;
    typedef typename Traits::FT FT;
    typedef typename Min_sphere_of_spheres_impl::Selector<FT>::Result Result;
    typedef typename Traits::Algorithm Algorithm;
    typedef typename std::vector<Sphere>::const_iterator
                     Support_sphere_iterator;
    typedef const Result *Center_coordinate_iterator;

  private: // internally used types:
    typedef typename Min_sphere_of_spheres_impl::Selector<FT>::Is_exact
                     Is_exact;
    typedef typename Traits::Coordinate_iterator C_it;
    typedef typename Traits::Use_square_roots Use_square_roots;
      
  private: // traits class:
    Traits t; // To allow the traits to not only vary at compile- but
              // also at runtime, we instantiate it here.
  
  public: // construction and destruction:
    inline Min_sphere_of_spheres_d(const Traits& traits = Traits());
  
    template<typename InputIterator>
    inline Min_sphere_of_spheres_d(InputIterator begin,InputIterator end,
                    const Traits& traits = Traits()) :
      t(traits), e(0), ss(t) {
      CGAL_assertion(ss.dim < 0);
      set(begin,end);
    }

  private: // internal helpers:
    inline void add(const Sphere& b);
  
    template<typename InputIterator>
    inline void add(InputIterator begin,InputIterator end) {
      while (begin != end) {
        add(*begin);
        ++begin;
      }
    }

  private: // forbid copy constructing and assignment (because our
           // pointers would be wrong then):
    Min_sphere_of_spheres_d(const Min_sphere_of_spheres_d&);
    Min_sphere_of_spheres_d& operator=(const Min_sphere_of_spheres_d&);
  
  public: // accessors:
    inline Support_sphere_iterator support_spheres_begin() const;
    inline Support_sphere_iterator support_spheres_end() const;
    inline int ambient_dimension() const;
    inline int dimension() const;
    inline const FT& discriminant() const;
    inline const Result& radius() const;
    inline Center_coordinate_iterator center_begin() const;
    inline Center_coordinate_iterator center_end() const;
    const Traits& traits() const;

  public: // predicates:
    inline bool is_empty() const;

  public: // modifiers:
    inline void clear() {
      ss.clear();
      ls.clear();
      l.clear();
      e = 0;
    }

    template<typename InputIterator>
    inline void set(InputIterator begin,InputIterator end) {
      clear();
      insert(begin,end);
    }

    template<typename InputIterator>
    inline void insert(InputIterator begin,InputIterator end) {
      add(begin,end);
      update();
    }

    inline void insert(const Sphere& s) {
      add(s);
      update();
    }
  
  public: // validity check:
    bool is_valid(bool verbose=false,int level=0) const;
  
  private:
    bool is_valid(bool verbose,const Tag_false is_exact) const;
    bool is_valid(bool verbose,const Tag_true is_exact) const;

    bool pivot(int B);

    void update();
    void update(LP_algorithm);
    void update(Farthest_first_heuristic);
  
    bool findFarthest(int from,int to,int& i,
		      const Tag_true use_sqrt,
		      const Tag_false is_exact);
    bool findFarthest(int from,int to,int& i,
		      const Tag_true use_sqrt,
		      const Tag_true is_exact);
    bool findFarthest(int from,int to,int& i,
		      const Tag_false use_sqrt,
		      const Tag_false is_exact);
    bool findFarthest(int from,int to,int& i,
		      const Tag_false use_sqrt,
		      const Tag_true is_exact);

  private:
    std::list<Sphere> ls;       // list of add()'ed balls
    std::vector<const Sphere*> l; // list of pointers to the balls in ls
    Min_sphere_of_spheres_impl::Support_set<Traits> ss; // current basis
    int e;                        // l[0..(e-1)] is basis
  };

  template<class Traits>
  Min_sphere_of_spheres_d<Traits>::
    Min_sphere_of_spheres_d(const Traits& traits) : t(traits), e(0), ss(t) {
    // ??? following will crash because ss.radius() is called for ss.dim==-1:
    // CGAL_assertion(Min_sphere_of_spheres_impl::isNeg(ss.radius(),ss.disc()));
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::add(const Sphere& b) {
    // allocate memory in case this is the first ball to add:
    if (ss.dim < 0) {
      CGAL_assertion_msg(l.size() == 0,"dimension not set but balls added");
      ss.reset(t.dimension(b));
    }

    // add ball:
    CGAL_precondition_msg(t.radius(b) >= 0,"sphere has negative radius");
    CGAL_precondition_msg(t.dimension(b)==ss.dim,"dimensions don't match");
    ls.push_back(b);
    l.push_back(&ls.back());
  }

  template<class Traits>
  const Traits& Min_sphere_of_spheres_d<Traits>::traits() const {
    return t;
  }

  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::is_empty() const {
    return Min_sphere_of_spheres_impl::isNeg(ss.radius(),ss.disc());
  }

  template<class Traits>
  int Min_sphere_of_spheres_d<Traits>::ambient_dimension() const {
    CGAL_assertion(ss.dim==-1 || !is_empty());
    return ss.dim;
  }

  template<class Traits>
  int Min_sphere_of_spheres_d<Traits>::dimension() const {
    CGAL_assertion(ss.dim==-1 || !is_empty());
    return is_empty()? ss.dim-1 : -1;
  }

  template<class Traits>
  const typename Min_sphere_of_spheres_d<Traits>::Result&
    Min_sphere_of_spheres_d<Traits>::radius() const {
    CGAL_assertion(!is_empty());
    return ss.radius();
  }

  template<class Traits>
  typename Min_sphere_of_spheres_d<Traits>::Center_coordinate_iterator
    Min_sphere_of_spheres_d<Traits>::center_begin() const {
    CGAL_assertion(!is_empty());
    return ss.begin();
  }

  template<class Traits>
  typename Min_sphere_of_spheres_d<Traits>::Center_coordinate_iterator
    Min_sphere_of_spheres_d<Traits>::center_end() const {
    CGAL_assertion(!is_empty());
    return ss.begin()+ss.dim;
  }

  template<class Traits>
  const typename Min_sphere_of_spheres_d<Traits>::FT&
    Min_sphere_of_spheres_d<Traits>::discriminant() const {
    CGAL_assertion(!is_empty());
    return ss.disc();
  }

  template<class Traits>
  inline void Min_sphere_of_spheres_d<Traits>::update() {
    // compute miniball:
    update(Algorithm());
  }

  template<class Traits>
  inline typename Min_sphere_of_spheres_d<Traits>::Support_sphere_iterator
    Min_sphere_of_spheres_d<Traits>::support_spheres_begin() const {
    return Support_sphere_iterator(l.begin());
  }

  template<class Traits>
  inline typename Min_sphere_of_spheres_d<Traits>::Support_sphere_iterator
    Min_sphere_of_spheres_d<Traits>::support_spheres_end() const {
    return Support_sphere_iterator(l.begin()+e);
  }

} // namespace CGAL

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Min_sphere_of_spheres_support_set.C>
#  include <CGAL/Min_sphere_of_spheres_pivot.C>
#  include <CGAL/Min_sphere_of_spheres_d.C>
#endif

#include <CGAL/Min_sphere_of_spheres_d_traits_d.h>

#endif // CGAL_MIN_SPHERE_OF_SPHERES_D_H
