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
// file          : include/CGAL/Min_sphere_of_spheres_d.C
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


#ifndef CGAL_MINIBALL_MINIBALL_C
#define CGAL_MINIBALL_MINIBALL_C

#include <numeric>
#include <CGAL/Min_sphere_of_spheres_d.h>

#ifdef CGAL_MINIBALL_DEBUG
#include <iostream>
#endif

#include <cstdlib>
#include <algorithm>
#include <CGAL/Min_sphere_of_spheres_d.h>

namespace CGAL_MINIBALL_NAMESPACE {

  // Returns a random number uniformly in {0,...,k-1}.
  inline int random(int k) {
    return std::rand() % k; // todo.  This is not good...
  }

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::update(LP_algorithm) {
    using namespace Min_sphere_of_spheres_d_impl;
    const int n = l.size();
    int i, k = n;
    do {
      CGAL_MINIBALL_ASSERT(k>=e && e>=0);
  
      // permute:
      for (int j=k-1; j>=e; --j) // todo. theory: is this necessary?
        std::swap(l[j],l[e+random(j+1-e)]);
  
      for (i=e; i<n; ++i)
        if (!ss.contains(t.begin(*l[i]),t.radius(*l[i]),Tol,Is_exact())) {
          pivot(i);
          k = i+1;
          break;
        }
    } while (i < n);
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::find_farthest(int from,int to,
    int& i,const Tag_true use_sqrt,const Tag_false is_exact) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // we will compute the excesses w.r.t. to the ball B with
    // center ss.begin() and radius radius:
    const FT radius = ss.radius() * Tol;
  
    // find ball with largest excess:
    FT maximum = radius;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const FT dist = inner_product_n<D>(ss.begin(),
        t.begin(*l[k]),0.0,std::plus<FT>(),
        Subtract_and_square<FT>());
  
      // compute excess:
      using std::sqrt;
      const FT ex = sqrt(dist)+t.radius(*l[k]);
  
      // compare with current maximum:
      if (ex > maximum) { // (*)
        maximum = ex;
        i = k;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return maximum > radius;
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::find_farthest(int from,int to,
    int& i,const Tag_true use_sqrt,const Tag_true is_exact) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // we will compute the excesses w.r.t. to the ball B with
    // center center and radius radius:
    Pair_to_double<FT> cast(ss.disc());
    const double radius = cast(ss.radius());
    double center[D];
    std::transform(ss.begin(),ss.begin()+D,center,cast);
  
    // find ball with largest excess:
    double maximum = radius;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const double dist = inner_product_n<D>(center,
         t.begin(*l[k]),0.0,std::plus<double>(),
         Subtract_and_square_to_double<FT>());
  
      // compute excess:
      using std::sqrt;
      const double ex = sqrt(dist) +
         CGAL_MINIBALL_NTS to_double(t.radius(*l[k]));
  
      // compare with current maximum:
      if (ex > maximum) { // (*)
        maximum = ex;
        i = k;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return maximum > radius &&
           !ss.contains(t.begin(*l[i]),t.radius(*l[i]),Tol,Is_exact());
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::find_farthest(int from,int to,
    int& i,const Tag_false use_sqrt,const Tag_false is_exact) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // we will compute the excesses w.r.t. to the ball with
    // center ss.begin() and radius radius:
    const FT radius = ss.radius() * Tol;
  
    // find ball with largest excess:
    bool found = false;
    FT max = radius, maxp = 0;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const FT dist = inner_product_n<D>(ss.begin(),
        t.begin(*l[k]),0.0,std::plus<FT>(),
        Subtract_and_square<FT>());
  
      if (compare(max,maxp,t.radius(*l[k]),dist)) {
        max   = t.radius(*l[k]);
        maxp  = dist;
        i     = k;
        found = true;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return found;
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::find_farthest(int from,int to,
    int& i,const Tag_false use_sqrt,const Tag_true is_exact) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // we will compute the excesses w.r.t. to the ball B with
    // center center and radius radius:
    Pair_to_double<FT> cast(ss.disc());
    const double radius = cast(ss.radius());
    double center[D];
    std::transform(ss.begin(),ss.begin()+D,center,cast);
  
    // find ball with largest excess:
    bool found = false;
    double max = radius, maxp = 0;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const double dist = inner_product_n<D>(center,
         t.begin(*l[k]),0.0,std::plus<double>(),
         Subtract_and_square_to_double<FT>());
  
      const double r = CGAL_MINIBALL_NTS to_double(t.radius(*l[k]));
      if (compare(max,maxp,r,dist)) {
        max   = r;
        maxp  = dist;
        i     = k;
        found = true;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return found &&
           !ss.contains(t.begin(*l[i]),t.radius(*l[i]),Tol,Is_exact());
  }
  
  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::update(Farthest_first_heuristic) {
    const int n = l.size();
    int i = e;
    CGAL_MINIBALL_ASSERT(e <= n);
  
    bool enclosing = (e == 0)? n == 0 :
      !find_farthest(e,n,i,Use_sqrt(),Is_exact());
  
    while (!enclosing && pivot(i)) {
      enclosing = !find_farthest(e,n,i,Use_sqrt(),Is_exact());
    }
  
    if (!is_approximate(Is_exact()))
      update(LP_algorithm());
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::
    is_valid(const Tag_false is_exact) {
    return true;
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::
    is_valid(const Tag_true is_exact) {
    using namespace Min_sphere_of_spheres_d_impl;
    using std::cerr;
    using std::endl;
  
    // check size of support set:
    if (e > l.size() || e > (D+1)) {
      cerr << "Min_sphere_of_spheres_d: support set too large." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    } else if (l.size() > 0 && e<=0) {
      cerr << "Min_sphere_of_spheres_d: support set too small." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
  
    // check case of no balls:
    if (l.size() <= 0 && !is_empty()) {
      cerr << "Min_sphere_of_spheres_d: miniball of {} non-empty." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
  
    // check that the miniball is enclosing:
    for (int i=0; i<l.size(); ++i)
      if (!ss.contains(t.begin(*l[i]),t.radius(*l[i]),Tol,Is_exact())) {
        cerr << "Min_sphere_of_spheres_d: miniball not enclosing." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // check that all support balls lie on the boundary:
    typedef Pair<FT> P;
    bool isSupporting = true;
    for (int i=0; i<e; ++i) {
      // check radii:
      const P rd = ss.radius()-t.radius(*l[i]);
      if (is_neg(rd,ss.disc()))
        isSupporting = false;
    
      // compute the (squared) distance from ss.begin() to l[i]'s center:
      const P dist = inner_product_n<D>(ss.begin(),
        t.begin(*l[i]),P(0,0),std::plus<P>(),
        Subtract_and_square_pair<FT>(ss.disc()));
    
      // compute the square of rd:
      const P sqrRd(sqr(rd.first)+sqr(rd.second)*ss.disc(),
                    FT(2)*rd.first*rd.second);
    
      // check containment:
      if (!is_zero(dist-sqrRd,ss.disc()))
        isSupporting = false;
    }
    if (!isSupporting) {
      cerr << "Min_sphere_of_spheres_d: support not on boundary." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
    
    // set up initial system's coefficient matrix:
    FT (*m)[D+1] = new FT[e][D+1];
    for (int j=0; j<e; ++j)
      copy_n<D>(t.begin(*l[j]),m[j]);
    for (int j=0; j<e; ++j)
      m[j][D] = FT(1);
    
    // set up initial system's right-hand-side:
    Pair<FT> rhs[D+1];
    copy_n<D>(ss.begin(),rhs);
    rhs[D] = FT(1);
    
    // perform Gaussian elimination:
    for (int j=0; j<e; ++j) {
      // check rank:
      if (m[j][j] == FT(0)) {
        // find row with non-zero entry in column j:
        int i = j;
        while (i<D+1 && m[j][i]==FT(0))
          ++i;
        if (i >= D+1) {
          cerr << "Min_sphere_of_spheres_d: supp. centers aff. dep." << endl
               << "Please contact the author <kf@iaeth.ch>." << endl;
          return false;
        }
    
        // exchange rows:
        for (int k=0; k<e; ++k)
          std::swap(m[k][j],m[k][i]);
        std::swap(rhs[j],rhs[i]);
      }
      CGAL_MINIBALL_ASSERT(m[j][j] != FT(0));
    
      // eliminate m[j][j+1..D] by subtracting a
      // multiple of row j from row i:
      for (int i=j+1; i<D+1; ++i) {
        // determine factor:
        const FT l = m[j][i]/m[j][j];
    
        // subtract row j times l from row i:
        for (int k=0; k<e; ++k)
          m[k][i] -= m[k][j]*l;
        rhs[i] -= rhs[j]*l;
      }
    }
    
    // check that we now have an upper triangular matrix:
    for (int j=0; j<e; ++j)
      for(int i=j+1; i<D+1; ++i)
        CGAL_MINIBALL_ASSERT(m[j][i] == FT(0));
    
    // check solvability:
    for (int i=e; i<D+1; ++i)
      if (!is_zero(rhs[i],ss.disc())) {
        cerr << "Min_sphere_of_spheres_d: center of the miniball" << endl
             << "not in the span of the support centers." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // compute coefficients by backsubstitution:
    Pair<FT> lambda[D+1];
    for (int i=e-1; i>=0; --i) {
      lambda[i] = rhs[i];
      for (int j=i+1; j<e; ++j)
        lambda[i] -= lambda[j]*m[j][i];
      lambda[i] = lambda[i]/m[i][i];
    }
    
    // check coefficients:
    for (int i=0; i<e; ++i)
      if (is_neg_or_zero(lambda[i],ss.disc())) {
        cerr << "Min_sphere_of_spheres_d: center of miniball not in" << endl
             << "interior of convex hull of support centers." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // tidy up:
    delete[] m;
    
    return true;
  }
  

} // namespace CGAL_MINIBALL_NAMESPACE

#endif // CGAL_MINIBALL_MINIBALL_CC
