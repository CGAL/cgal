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
// file          : include/CGAL/Min_sphere_of_spheres_support_set.h
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

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_H
#define CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_H

#include <functional>
#include <cmath>

#include <CGAL/basic.h>
#include <CGAL/Min_sphere_of_spheres_pair.h>

namespace CGAL {

  // another namespace to "hide" implementation details:
  namespace Min_sphere_of_spheres_impl {

    // Constants used for approximate computation (i.e. when FT is double):
    const double Min = 1.0e-120;
    const double Eps = 1.0e-16;
    const double Tol = 1.0+Eps;
    
    // The following is used to construct at compile-time the type of
    // the radius and the coordinates of the miniball.  For
    // approximate computation this is simply the floating-point type
    // FT, for exact computation it is Pair<FT>.
    
    template<typename FT>
    struct Selector {
      typedef Pair<FT> Result;
      typedef Tag_true Is_exact;
    };
  
    template<>
    struct Selector<double> {
      typedef double Result;
      typedef Tag_false Is_exact;
    };

    inline bool is_approximate(const Tag_false is_exact) {
      return true;
    }

    inline bool is_approximate(const Tag_true is_exact) {
      return false;
    }

    // For the routine Support_set<Traits>::contains() we need the
    // following functors (in order to formulate the code using
    // std::inner_product()):
    
    template<typename FT>
    struct Subtract_and_square {
      inline FT operator()(const FT x,const FT y) const {
	return square(x-y);
      }
    };
    
    template<typename FT>
    class Subtract_and_square_pair {
    private:
      typedef Min_sphere_of_spheres_impl::Pair<FT> Pair;
      const FT& d;
      
    public:
      Subtract_and_square_pair(const FT& d) : d(d) {}
      
      inline Pair operator()(const Pair& x,const FT& y) const {
	const FT t = x.first - y;
	return Pair(CGAL_NTS square(t) + CGAL_NTS square(x.second)*d,
		    // Note: The following line is what I would like
		    // to write.  However, gcc 2.95.2 fails to
		    // implicitly convert 2 to FT, at least for FT
		    // being Gmpq (which doesn't have a specialized
		    // version of operator* for pairs (int,Gmpq).
		    // 2*t*x.second);
		    static_cast<FT>(2)*t*x.second);
      }
    };
    
    template<typename FT>
    class Pair_to_double {
    private:
      typedef Min_sphere_of_spheres_impl::Pair<FT> Pair;
      const double root;
      
    public:
      Pair_to_double(const FT& disc) :
	root(std::sqrt(CGAL_NTS to_double(disc))) {}
      
      inline double operator()(const Pair& p) const {
	return CGAL_NTS to_double(p.first) +
               CGAL_NTS to_double(p.second) * root;
      }
    };
  
    template<typename FT>
    struct Subtract_and_square_to_double {
      inline double operator()(const double x,const FT& y) const {
	return CGAL_NTS square(x - CGAL_NTS to_double(y));
      }
    };

    template<class Traits>
    class Support_set {
    public: // public types:
      typedef typename Traits::Sphere Sphere;
      typedef typename Traits::FT FT;
      typedef typename Selector<FT>::Result Result;
      typedef typename Traits::Algorithm Algorithm;
      typedef typename Selector<FT>::Is_exact Is_exact;

    private: // internally used types:
      typedef typename Traits::Coordinate_iterator C_it;
      
    public: // constructor:
      inline Support_set(Traits& traits);
      ~Support_set();
      
    public: // access:
      inline const Result& radius() const;
      inline const Result *begin() const;
      inline const FT& disc() const;
      
    public: // containment test:
      template<typename InputIterator>
      bool contains(InputIterator c,const FT& r,
		    const double tol,
		    const Tag_false is_exact) const {
	// scale ball:
	const FT r1 = sol[k]*tol;
	
	// check radii:
	if (r > r1)
	  return false;
	
	// compute the (squared) distance from center[k] to c:
	const FT dist = std::inner_product(center[k],center[k]+dim,c,
	  FT(0),std::plus<FT>(),Subtract_and_square<FT>());
    
	// check containment:
	return dist <= CGAL_NTS square(r1-r);
      }
    
      template<typename InputIterator>
      bool contains(InputIterator c,const FT& r,
		    const double,const Tag_true is_exact) const {
	typedef Pair<FT> P;
	
	// check radii:
	const P rd = sol[k]-r;
	if (isNeg(rd,discrim))
	  return false;
	
	// compute the (squared) distance from center[k] to c:
	P dist = std::inner_product(center[k],center[k]+dim,
          c,P(0,0),std::plus<P>(),Subtract_and_square_pair<FT>(discrim));
	Min_sphere_of_spheres_impl::normalize(dist);
    
	// compute the square of rd:
	const P sqrRd(CGAL_NTS square(rd.first) +
		      CGAL_NTS square(rd.second)*discrim,
		      // Note: The following line is what I would like
		      // to write.  However, gcc 2.95.2 fails to
		      // implicitly convert 2 to FT, at least for FT
		      // being Gmpq (which doesn't have a specialized
		      // version of operator* for pairs (int,Gmpq).
		      // 2*rd.first*rd.second);
		      static_cast<FT>(2)*rd.first*rd.second);
	
	// check containment:
	return isNegOrZero(dist-sqrRd,discrim);
      }

    public: // initialization:
      void reset(int dim);
      
    public: // modification:
      void clear();
      bool push(const Sphere& ball);
      inline void pop();
      bool isValid();
  
    private: // utility:
      void allocate(int dim);
      void deallocate();
      inline int findRadii(const Tag_false is_exact);
      inline int findRadii(const Tag_true is_exact);
      inline void findCenters(int nr);
      inline bool selectBall(int nr);
      inline bool isSupporting(int i) const;
      
    private: // traits class:
      Traits& t;

    public: // member fields:
      int dim;                  // dimension of the balls (or -1 if unknown)

    private: // member fields:
      int m;                    // FT of pushed balls
      const Sphere **b;         // pointers to pushed balls
      
      // variables of the device:
      FT **u;
      FT **d;
      FT **e;
      FT **f;
      FT *alpha;
      FT *delta;
      FT *eps;
      FT *phi;
      FT *sigma;
      FT *chi;
      FT *psi;
      FT *omega;
      FT **tau;
      
      FT discrim;          // discriminant of the quadratic
      Result sol[2];       // solutions of the quadric (sol[0] <= sol[1])
      Result **center;     // two centers corresponding to the above radii
      Result **beta;       // beta[i] is equal to $\beta_{i+1}$
      int k;

    public:
      Result *temp;        // used temporarily in isSupporting()
      double *float_temp;  // used temporarily in findFarthest()
    };
    
    template<class Traits>
    Support_set<Traits>::Support_set(Traits& traits) : t(traits), dim(-1) {}

    template<class Traits>
    const typename Support_set<Traits>::Result&
      Support_set<Traits>::radius() const {
      CGAL_assertion(dim != -1);
      return sol[k];
    }

    template<class Traits>
    const typename Support_set<Traits>::Result
      *Support_set<Traits>::begin() const {
      CGAL_assertion(dim != -1);
      return center[k];
    }
    
    template<class Traits>
    const typename Support_set<Traits>::FT&
      Support_set<Traits>::disc() const {
      CGAL_assertion(dim != -1);
      return discrim;
    }
    
    template<class Traits>
    void Support_set<Traits>::pop() {
      CGAL_assertion(dim != -1);
      CGAL_assertion(m > 0);
      --m;
    }
  
    template<class Traits>
    int Support_set<Traits>::findRadii(const Tag_false is_exact) {
      CGAL_assertion(dim != -1);

      // find discriminant:
      const double disc = CGAL_NTS square(psi[m-1]) - 4*chi[m-1]*omega[m-1];
      if (disc < Min)
	return 0;
      
      // handle degenerate quadratic:
      if (std::abs(chi[m-1]) < Min) {
	if (std::abs(psi[m-1]) < Min)
	  return 0;
	sol[0] = -omega[m-1]/psi[m-1];
	return 1;
      }
      
      // calculate the two solutions:
      double sd = std::sqrt(disc);
      if (psi[m-1] > 0)
	sd = -sd;
      sol[0] = 2*omega[m-1] / (sd-psi[m-1]);
      sol[1] = (sd-psi[m-1]) / (2*chi[m-1]);
      
      // order the solutions (*):
      if (sol[0] > sol[1])
	std::swap(sol[0],sol[1]);
      
      // eliminate negative solutions:
      if (sol[0] < 0) {
	if (sol[1] < 0)
	  return 0;
	sol[0] = sol[1];
	return 1;
      }
      return 2;
    }
    
    template<class Traits>
    int Support_set<Traits>::findRadii(const Tag_true is_exact) {
      CGAL_assertion(dim != -1);

      // find discriminant:
      // Note: The following line is what I would like
      // to write.  However, gcc 2.95.2 fails to
      // implicitly convert 2 to FT, at least for FT
      // being Gmpq (which doesn't have a specialized
      // version of operator* for pairs (int,Gmpq).
      // discrim = CGAL_NTS square(psi[m-1]) - 4*chi[m-1]*omega[m-1];
      discrim = CGAL_NTS square(psi[m-1]) - FT(4)*chi[m-1]*omega[m-1];
      if (discrim <= 0)
	return 0;
      
      // handle degenerate quadratic:
      if (chi[m-1] == 0) {
	if (psi[m-1] == 0)
	  return 0;
	sol[0] = Pair<FT>(-omega[m-1]/psi[m-1],0);
	return 1;
      }
      
      // calculate the two solutions:
      // Note: The following two lines is what I would like to write.
      // However, gcc 2.95.2 fails to implicitly convert 2 to FT, at
      // least for FT being Gmpq (which doesn't have a specialized
      // version of operator* for pairs (int,Gmpq).
      // sol[0] = Pair<FT>(-psi[m-1]/(2*chi[m-1]),-1/(2*chi[m-1]));
      // sol[1] = Pair<FT>(-psi[m-1]/(2*chi[m-1]),+1/(2*chi[m-1]));
      sol[0] = Pair<FT>(-psi[m-1]/(FT(2)*chi[m-1]),FT(-1)/(FT(2)*chi[m-1]));
      sol[1] = Pair<FT>(-psi[m-1]/(FT(2)*chi[m-1]),FT(+1)/(FT(2)*chi[m-1]));
      
      // order the solutions (*):
      if (chi[m-1] < 0)
	std::swap(sol[0],sol[1]);
      CGAL_assertion(isNegOrZero(sol[0]-sol[1],discrim));
      
      // eliminate negative solutions:
      if (isNeg(sol[0],discrim)) {
	if (isNeg(sol[1],discrim))
	  return 0;
	sol[0] = sol[1];
	return 1;
      }
      return 2;
    }
  
    template<class Traits>
    void Support_set<Traits>::findCenters(int) {
      CGAL_assertion(dim != -1);

      // initialize both centers to the center of the first pushed ball:
      C_it c = t.center_coordinates_begin(*b[0]);
      for(int j=0; j<dim; ++j, ++c)
	center[0][j] = center[1][j] = *c;
      
      // compute centers:
      for(int i=1; i<m; ++i) {
	beta[0][i-1] = (delta[i]+eps[i]+sol[0]*phi[i])/alpha[i];
	beta[1][i-1] = (delta[i]+eps[i]+sol[1]*phi[i])/alpha[i];
	Min_sphere_of_spheres_impl::normalize(beta[0][i-1]);
	Min_sphere_of_spheres_impl::normalize(beta[1][i-1]);
	for (int j=0; j<dim; ++j) {
	  center[0][j] += beta[0][i-1]*u[i][j];
	  center[1][j] += beta[1][i-1]*u[i][j];
	  Min_sphere_of_spheres_impl::normalize(center[0][j]);
	  Min_sphere_of_spheres_impl::normalize(center[1][j]);
	}
      }
    }
    
    template<typename FT>
    inline bool isPoint(const FT& radius,const double tol,
			const FT& sol,const Tag_false is_exact) {
      return radius > tol*sol;
    }
    
    template<typename FT>
    inline bool isPoint(const FT& radius,const double,
			const Pair<FT>&,const Tag_true is_exact) {
      return radius > 0;
    }
    
    template<class Traits>
    bool Support_set<Traits>::selectBall(int nr) {
      CGAL_assertion(dim != -1);
      CGAL_assertion(nr==1 || nr==2);
  
      // take the first (ie. smallest) ball which is enclosing:
      for (k=0; k<nr; ++k) {
	bool isEnclosing = true;
	for (int i=0; i<m; ++i)
	  if (isPoint(t.radius(*b[i]),Eps,sol[k],Is_exact()))
	    if (isNeg(sol[k]-t.radius(*b[i]),discrim) ||
		!contains(t.center_coordinates_begin(*b[i]),
			  FT(0),1,Is_exact())) {
	      isEnclosing = false;
	      break;
	    }
	
	if (isEnclosing)
	  return isSupporting(k);
      }
      
      // signal ``no'' if no ball was enclosing:
      return false;
    }
  
    template<class Traits>
    bool Support_set<Traits>::isSupporting(int k) const {
      CGAL_assertion(dim != -1);

      Result ming(0);
      Result gamma0(1);
      
      for (int i=m-1; i>0; --i) {
	temp[i-1] = beta[k][i-1];
	for (int j=i+1; j<m; ++j)
	  temp[i-1] -= temp[j-1]*tau[i][j];
	gamma0 -= temp[i-1];
	Min_sphere_of_spheres_impl::normalize(temp[i-1]);
	Min_sphere_of_spheres_impl::normalize(gamma0);
	if (isNeg(temp[i-1]-ming,discrim))
	  ming = temp[i-1];
      }
      if (isNeg(gamma0-ming,discrim))
	ming = gamma0;

      return !isNeg(ming,discrim);
    }
  
  } // namespace Mini_sphere_of_spheres_impl

} // namespace CGAL

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Min_sphere_of_spheres_support_set.C>
#endif

#endif // CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_H
