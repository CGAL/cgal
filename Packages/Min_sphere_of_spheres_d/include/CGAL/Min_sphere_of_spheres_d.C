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
// file          : include/CGAL/Min_sphere_of_spheres.C
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

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_D_C
#define CGAL_MIN_SPHERE_OF_SPHERES_D_C

#include <numeric>
#include <cstdlib>
#include <algorithm>

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Random.h>

namespace CGAL {

  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::update(const LP_algorithm) {
    CGAL::Random random;
    const int n = l.size();
    int i, k = n;
    do {
      CGAL_assertion(k>=e && e>=0);
  
      // permute:
      for (int j=k-1; j>=e; --j) // todo. theory: is this necessary?
        std::swap(l[j],l[random.get_int(e,j+1)]);
  
      for (i=e; i<n; ++i) // todo. not optimal for updates
        if (!ss.contains(t.center_coordinates_begin(*l[i]),t.radius(*l[i]),
			 Min_sphere_of_spheres_impl::Tol,Is_exact())) {
          pivot(i);
          k = i+1;
          break;
        }
    } while (i < n);
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::findFarthest(int from,int to,int& i,
    const Tag_true use_sqrt,const Tag_false is_exact) {
    // we will compute the excess's wrt. to the ball B with
    // center ss.begin() and radius radius:
    const FT radius = ss.radius() * Min_sphere_of_spheres_impl::Tol;
  
    // find ball with largest excess:
    FT maximum = radius;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const FT dist = std::inner_product(ss.begin(),ss.begin()+ss.dim,
        t.center_coordinates_begin(*l[k]),0.0,std::plus<FT>(),
        Min_sphere_of_spheres_impl::Subtract_and_square<FT>());

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
  bool Min_sphere_of_spheres_d<Traits>::findFarthest(int from,int to,int& i,
    const Tag_true use_sqrt,const Tag_true is_exact) {
    // we will compute the excess's wrt. to the ball B with
    // center ss.float_temp and radius radius:
    Min_sphere_of_spheres_impl::Pair_to_double<FT> cast(ss.disc());
    const double radius = cast(ss.radius());
    std::transform(ss.begin(),ss.begin()+ss.dim,ss.float_temp,cast);
  
    // find ball with largest excess:
    double maximum = radius;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const double dist = std::inner_product(ss.float_temp,
         ss.float_temp+ss.dim,t.center_coordinates_begin(*l[k]),
         0.0,std::plus<double>(),
         Min_sphere_of_spheres_impl::Subtract_and_square_to_double<FT>());
  
      // compute excess:
      using std::sqrt;
      const double ex = sqrt(dist)+to_double(t.radius(*l[k]));
  
      // compare with current maximum:
      if (ex > maximum) { // (*)
        maximum = ex;
        i = k;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return maximum > radius &&
           !ss.contains(t.center_coordinates_begin(*l[i]),t.radius(*l[i]),
			Min_sphere_of_spheres_impl::Tol,Is_exact());
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::findFarthest(int from,int to,int& i,
    const Tag_false use_sqrt,const Tag_false is_exact) {
    // we will compute the excess's wrt. to the ball with
    // center ss.begin() and radius radius:
    const FT radius = ss.radius() * Min_sphere_of_spheres_impl::Tol;
  
    // find ball with largest excess:
    bool found = false;
    FT max = radius, maxp = 0;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const FT dist = std::inner_product(ss.begin(),ss.begin()+ss.dim,
        t.center_coordinates_begin(*l[k]),0.0,std::plus<double>(),
        Min_sphere_of_spheres_impl::Subtract_and_square<FT>());
  
      if (Min_sphere_of_spheres_impl::compare(max,maxp,t.radius(*l[k]),dist)) {
        max  = t.radius(*l[k]);
        maxp = dist;
        i = k;
        found = true;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return found;
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::findFarthest(int from,int to,int& i,
    const Tag_false use_sqrt,const Tag_true is_exact) {
    // we will compute the excess's wrt. to the ball B with
    // center ss.float_temp and radius radius:
    Min_sphere_of_spheres_impl::Pair_to_double<FT> cast(ss.disc());
    const double radius = cast(ss.radius());
    std::transform(ss.begin(),ss.begin()+ss.dim,ss.float_temp,cast);
  
    // find ball with largest excess:
    bool found = false;
    double max = radius, maxp = 0;
    for (int k=from; k<to; ++k) {
      // compute the (squared) distance from c1 to c2:
      const double dist = std::inner_product(ss.float_temp,
         ss.float_temp+ss.dim,t.center_coordinates_begin(*l[k]),
         0.0,std::plus<double>(),
         Min_sphere_of_spheres_impl::Subtract_and_square_to_double<FT>());
  
      const double r = CGAL_NTS to_double(t.radius(*l[k]));
      if (Min_sphere_of_spheres_impl::compare(max,maxp,r,dist)) {
        max  = r;
        maxp = dist;
        i = k;
        found = true;
      }
    }
  
    // return whether B doesn't contain the ball l[i]:
    return found &&
           !ss.contains(t.center_coordinates_begin(*l[i]),t.radius(*l[i]),
			Min_sphere_of_spheres_impl::Tol,Is_exact());
  }
  
  template<class Traits>
  void Min_sphere_of_spheres_d<Traits>::update(Farthest_first_heuristic) {
    const int n = l.size();
    int i = e;
    CGAL_assertion(e <= n);
  
    bool enclosing = (e == 0)? n == 0 :
      !findFarthest(e,n,i,Use_square_roots(),Is_exact());
  
    while (!enclosing && pivot(i))
      enclosing = !findFarthest(e,n,i,Use_square_roots(),Is_exact());

    if (!Min_sphere_of_spheres_impl::is_approximate(Is_exact()))
      update(LP_algorithm());
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::is_valid(bool verbose,int) const {
    return is_valid(verbose,Is_exact());
  }

  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::
    is_valid(bool verbose,const Tag_false is_exact) const {
    // ???
    // CGAL_warning_msg(false,
    // "is_valid() not supported for inexact number types");
    return true;
  }

  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::
    is_valid(bool verbose,const Tag_true is_exact) const {
    typedef Min_sphere_of_spheres_impl::Pair<FT> P;
    Verbose_ostream verr(verbose);
    using std::endl;
  
    // check size of support set:
    if (e > l.size() || e > (ss.dim+1)) {
      verr << "Package Miniball: support set too large." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    } else if (l.size() > 0 && e<=0) {
      verr << "Package Miniball: support set too small." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
  
    // check case of no balls:
    if (l.size() <= 0 && !is_empty()) {
      verr << "Package Miniball: miniball of {} is non-empty." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
  
    // check that the miniball is enclosing:
    for (int i=0; i<l.size(); ++i)
      if (!ss.contains(t.center_coordinates_begin(*l[i]),t.radius(*l[i]),
		       Min_sphere_of_spheres_impl::Tol,Is_exact())) {
        verr << "Package Miniball: miniball not enclosing." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // check that all support balls lie on the boundary:
    bool isSupporting = true;
    for (int i=0; i<e; ++i) {
      // check radii:
      const P rd = ss.radius()-t.radius(*l[i]);
      if (isNeg(rd,ss.disc()))
        isSupporting = false;
    
      // compute the (squared) distance from ss.begin() to l[i]'s center:
      P dist = std::inner_product(ss.begin(),ss.begin()+ss.dim,
        t.center_coordinates_begin(*l[i]),P(0,0),std::plus<P>(),
        Min_sphere_of_spheres_impl::Subtract_and_square_pair<FT>(ss.disc()));
      Min_sphere_of_spheres_impl::normalize(dist); 
    
      // compute the square of rd:
      P sqrRd(CGAL_NTS square(rd.first) + CGAL_NTS square(rd.second)*ss.disc(),
	      // Note: The following line is what I would like
	      // to write.  However, gcc 2.95.2 fails to
	      // implicitly convert 2 to FT, at least for FT
	      // being Gmpq (which doesn't have a specialized
	      // version of operator* for pairs (int,Gmpq).
	      // 2*rd.first*rd.second);
	      FT(2)*rd.first*rd.second);
    
      // check containment:
      if (!Min_sphere_of_spheres_impl::isZero(dist-sqrRd,ss.disc()))
        isSupporting = false;
    }
    if (!isSupporting) {
      verr << "Package Miniball: support ball not on boundary." << endl
           << "Please contact the author <kf@iaeth.ch>." << endl;
      return false;
    }
    
    // set up initial system's coefficient matrix:
    FT **m = new FT*[e];
    for (int j=0; j<e; ++j)
      m[j] = new FT[ss.dim+1];
    for (int j=0; j<e; ++j)
      std::copy(t.center_coordinates_begin(*l[j]),
		t.center_coordinates_begin(*l[j])+ss.dim,m[j]);
    for (int j=0; j<e; ++j)
      m[j][ss.dim] = FT(1);
    
    // set up initial system's right-hand-side:
    P *rhs = new P[ss.dim+1];
    std::copy(ss.begin(),ss.begin()+ss.dim,rhs);
    rhs[ss.dim] = FT(1);
    
    // perform Gaussian elimination:
    for (int j=0; j<e; ++j) {
      // check rank:
      if (m[j][j] == 0) {
        // find row with non-zero entry in column j:
        int i = j;
        while (i<ss.dim+1 && m[j][i]==0)
          ++i;
        if (i >= ss.dim+1) {
          verr << "Package Miniball: supp. centers not aff. indep." << endl
               << "Please contact the author <kf@iaeth.ch>." << endl;
          return false;
        }
    
        // exchange rows:
        for (int k=0; k<e; ++k)
          std::swap(m[k][j],m[k][i]);
        std::swap(rhs[j],rhs[i]);
      }
      CGAL_assertion(m[j][j] != 0);
    
      // eliminate m[j][j+1..ss.dim] by subtracting a
      // multiple of row j from row i:
      for (int i=j+1; i<ss.dim+1; ++i) {
        // determine factor:
	FT l = m[j][i]/m[j][j];
	Min_sphere_of_spheres_impl::normalize(l);
    
        // subtract row j times l from row i:
        for (int k=0; k<e; ++k) {
          m[k][i] -= m[k][j]*l;
          Min_sphere_of_spheres_impl::normalize(m[k][i]);
        }
        rhs[i] -= rhs[j]*l;
	Min_sphere_of_spheres_impl::normalize(rhs[i]);
      }
    }
    
    // check that we now have an upper triangular matrix:
    for (int j=0; j<e; ++j)
      for(int i=j+1; i<ss.dim+1; ++i)
        CGAL_assertion(m[j][i] == 0);
    
    // check solvability:
    for (int i=e; i<ss.dim+1; ++i)
      if (!Min_sphere_of_spheres_impl::isZero(rhs[i],ss.disc())) {
        verr << "Package Miniball: center of the miniball not in" << endl
             << "the span of the support centers." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // compute coefficients by backsubstitution:
    P *lambda = new P[ss.dim+1];
    for (int i=e-1; i>=0; --i) {
      lambda[i] = rhs[i];
      for (int j=i+1; j<e; ++j) {
        lambda[i] -= lambda[j]*m[j][i];
	Min_sphere_of_spheres_impl::normalize(lambda[i]);
      }
      lambda[i] = lambda[i]/m[i][i];
    }
    
    // check coefficients:
    for (int i=0; i<e; ++i)
      if (Min_sphere_of_spheres_impl::isNegOrZero(lambda[i],ss.disc())) {
        verr << "Package Miniball: center of miniball not in" << endl
             << "interior of convex hull of support centers." << endl
             << "Please contact the author <kf@iaeth.ch>." << endl;
        return false;
      }
    
    // tidy up:
    delete[] rhs;
    delete[] lambda;
    for (int j=0; j<e; ++j)
      delete[] m[j];
    delete[] m;

    return true;
  }
  
} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_D_C
