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
// file          : include/CGAL/Min_sphere_of_spheres_support_set.C
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

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_C
#define CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_C

#include <CGAL/Min_sphere_of_spheres_support_set.h>

namespace CGAL {

  // another namespace to "hide" implementation details:
  namespace Min_sphere_of_spheres_impl {
    
    template<class Traits>
    Support_set<Traits>::~Support_set() {
      if (dim >= 0)
	deallocate();
    }
    
    template<typename FT>
    FT **allocate_matrix(const int m,const int n) {
      FT **p = new FT*[m];
      for (int i=0; i<m; ++i)
	p[i] = new FT[n];
      return p;
    }

    template<typename FT>
    void deallocate_matrix(FT** p,const int m,const int n) {
      for (int i=0; i<m; ++i)
	delete[] p[i];
      delete[] p;
    }

    template<class Traits>
    void Support_set<Traits>::allocate(const int dim) {
      CGAL_assertion(this->dim==-1 && dim>=0);
      b = new const Sphere*[dim+1];

      u = allocate_matrix<FT>(dim+1,dim);
      d = allocate_matrix<FT>(dim+1,dim);
      e = allocate_matrix<FT>(dim+1,dim);
      f = allocate_matrix<FT>(dim+1,dim);
      alpha = new FT[dim+1];
      delta = new FT[dim+1];
      eps   = new FT[dim+1];
      phi   = new FT[dim+1];
      sigma = new FT[dim+1];
      chi   = new FT[dim+1];
      psi   = new FT[dim+1];
      omega = new FT[dim+1];
      tau = allocate_matrix<FT>(dim,dim+1);
      
      center = allocate_matrix<Result>(2,dim);
      beta   = allocate_matrix<Result>(2,dim);
      temp   = new Result[dim+1];
      float_temp = new double[dim+1];
      this->dim = dim;
    }

    template<class Traits>
    void Support_set<Traits>::deallocate() {
      CGAL_assertion(dim != -1);

      delete[] b;

      deallocate_matrix(u,dim+1,dim);
      deallocate_matrix(d,dim+1,dim);
      deallocate_matrix(e,dim+1,dim);
      deallocate_matrix(f,dim+1,dim);
      delete[] alpha;
      delete[] delta;
      delete[] eps;
      delete[] phi;
      delete[] sigma;
      delete[] chi;
      delete[] psi;
      delete[] omega;
      delete[] tau;
      
      deallocate_matrix(center,2,dim);
      deallocate_matrix(beta  ,2,dim);
      delete[] temp;
      dim = -1;
    }

    template<class Traits>
    void Support_set<Traits>::reset(int dim) {
      // (re)allocate resources:
      CGAL_assertion(dim != -1);
      if (this->dim != -1)
	deallocate();
      allocate(dim);

      clear();
    }

    template<class Traits>
    void Support_set<Traits>::clear() {
      m = 0;
      k = 0;
      sol[0] = FT(-1);
    }

    template<typename FT>
    inline bool reject(const FT& alpha,const double tol,
		       const Tag_false is_exact) {
      return CGAL_NTS abs(alpha) < tol;
    }

    template<typename FT>
    inline bool reject(const FT& alpha,const double,
		       const Tag_true is_exact) {
      return alpha == 0;
    }

    template<class Traits>
    bool Support_set<Traits>::push(const Sphere& ball) {
      CGAL_assertion(dim != -1);

      // reject if the current basis has already dim+1 balls:
      if (m > dim)
	return false;
      
      // keep a pointer to the ball:
      b[m] = &ball;

      // handle first pushed ball:
      if (m == 0) {
	for (int j=0; j<dim; ++j)
	  d[0][j] = e[0][j] = f[0][j] = 0;
	sigma[0] = 0;
	chi[0]   = -2;
	// Note: The following line is what I would like
	// to write.  However, gcc 2.95.2 fails to
	// implicitly convert 2 to FT, at least for FT
	// being Gmpq (which doesn't have a specialized
	// version of operator* for pairs (int,Gmpq).
	// psi[0]   =  4 * t.radius(ball);
	psi[0]   =  static_cast<FT>(4) * t.radius(ball);
	// Note: The following line is what I would like
	// to write.  However, gcc 2.95.2 fails to
	// implicitly convert 2 to FT, at least for FT
	// being Gmpq (which doesn't have a specialized
	// version of operator* for pairs (int,Gmpq).
	// omega[0] = -2 * square(t.radius(ball));
	omega[0] = static_cast<FT>(-2) * square(t.radius(ball));
      }
      
      // handle second, third, etc. push:
      else {
	// calculate $C_m$, storing it (temporarily) in u[m]:
	C_it c = t.center_coordinates_begin(*b[m]),
	  c_0  = t.center_coordinates_begin(*b[0]);
	for (int j=0; j<dim; ++j) {
	  u[m][j] = *c - *c_0;
	  ++c; ++c_0;
	}
	
	// compute $\tau_{im}$ for $1<=i<m$:
	for (int i=1; i<m; ++i) {
	  tau[i][m] = 0;
	  for (int j=0; j<dim; ++j)
	    tau[i][m] += u[i][j]*u[m][j];
	  // Note: The following line is what I would like
	  // to write.  However, gcc 2.95.2 fails to
	  // implicitly convert 2 to FT, at least for FT
	  // being Gmpq (which doesn't have a specialized
	  // version of operator* for pairs (int,Gmpq).
	  // tau[i][m] *= 2/alpha[i];
	  tau[i][m] *= static_cast<FT>(2)/alpha[i];
	  Min_sphere_of_spheres_impl::normalize(tau[i][m]);
	}
	
	// calculate delta[m], eps[m] and phi[m] (by definition): (We might
	// compute these values in vain here (if alpha[m] gets too small
	// below).  Nonetheless we do the computation here because we have
	// $C_m$ at hand (namely in u[m]).)
	const FT t1 = t.radius(*b[0]) - t.radius(*b[m]),
                 t2 = t.radius(*b[0]) + t.radius(*b[m]);
	phi[m]   = eps[m] = 0;
	delta[m] = -sigma[m-1];
	for (int j=0; j<dim; ++j) {
	  eps[m]   -= u[m][j]*e[m-1][j];
	  phi[m]   -= u[m][j]*f[m-1][j];
	  delta[m] += square(u[m][j]-d[m-1][j]);
	}
	// Note: The following line is what I would like
	// to write.  However, gcc 2.95.2 fails to
	// implicitly convert 2 to FT, at least for FT
	// being Gmpq (which doesn't have a specialized
	// version of operator* for pairs (int,Gmpq).
	// phi[m] = 2*(phi[m] - t1);
	phi[m] = static_cast<FT>(2)*(phi[m] - t1);
	// Note: The following line is what I would like
	// to write.  However, gcc 2.95.2 fails to
	// implicitly convert 2 to FT, at least for FT
	// being Gmpq (which doesn't have a specialized
	// version of operator* for pairs (int,Gmpq).
	// eps[m] = t1*t2+2*eps[m];
	eps[m] = t1*t2+static_cast<FT>(2)*eps[m];
	Min_sphere_of_spheres_impl::normalize(eps[m]);
	Min_sphere_of_spheres_impl::normalize(phi[m]);
	Min_sphere_of_spheres_impl::normalize(delta[m]);
	  
	// fix u[m] to be $C_m-\q{C}_m$:
	// (This is only necessary for m>1 because $\q{C}_1=0$.)
	for (int i=1; i<m; ++i)
	  for (int j=0; j<dim; ++j) {
	    u[m][j] -= tau[i][m]*u[i][j];
	    Min_sphere_of_spheres_impl::normalize(u[m][j]);
	  }
	
	// calculate alpha[m]:
	alpha[m] = 0;
	for (int j=0; j<dim; ++j)
	  alpha[m] += square(u[m][j]);
	alpha[m] *= 2;
	Min_sphere_of_spheres_impl::normalize(alpha[m]);
	
	// reject push if alpha[m] is to small:
	if (reject(alpha[m],Eps,Is_exact()))
	  return false;
	
	// calculate d[m], e[m] and f[m]:
	const FT da = delta[m]/alpha[m],
	  ea = eps[m]/alpha[m],
	  fa = phi[m]/alpha[m];
	for (int j=0; j<dim; ++j) {
	  d[m][j] = d[m-1][j] + da*u[m][j];
	  e[m][j] = e[m-1][j] + ea*u[m][j];
	  f[m][j] = f[m-1][j] + fa*u[m][j];
	  Min_sphere_of_spheres_impl::normalize(d[m][j]);
	  Min_sphere_of_spheres_impl::normalize(e[m][j]);
	  Min_sphere_of_spheres_impl::normalize(f[m][j]);
	}
	
	// compute sigma[m], chi[m], psi[m] and omega[m]:
	const FT de = delta[m]+eps[m];
	sigma[m] = sigma[m-1] + delta[m]*da/2;
	chi[m]   = chi[m-1]   + phi[m]*fa;
        // Note: The following line is what I would like
        // to write.  However, gcc 2.95.2 fails to
        // implicitly convert 2 to FT, at least for FT
        // being Gmpq (which doesn't have a specialized
        // version of operator* for pairs (int,Gmpq).
	// psi[m]   = psi[m-1]   + 2*de*fa;
	psi[m]   = psi[m-1]   + static_cast<FT>(2)*de*fa;
	omega[m] = omega[m-1] + de*de/alpha[m];
	Min_sphere_of_spheres_impl::normalize(sigma[m]);
	Min_sphere_of_spheres_impl::normalize(chi[m]);
	Min_sphere_of_spheres_impl::normalize(psi[m]);
	Min_sphere_of_spheres_impl::normalize(omega[m]);
      }
      
      ++m;
      return true;
    }
    
    template<class Traits>
    bool Support_set<Traits>::isValid() {
      CGAL_assertion(dim != -1);
      CGAL_assertion(m > 0);            // we want at least one pushed ball
  
      if (m == 1) {
	std::copy(t.center_coordinates_begin(*b[0]),
		  t.center_coordinates_begin(*b[0])+dim,center[0]);
	sol[0] = t.radius(*b[0]);
	CGAL_assertion(k == 0);
	
      } else {
	// solve the quadratic (which yields two "candidate" balls):
	const int nr = findRadii(Is_exact());
	if (nr == 0)
	  return false;
	
	// compute centers:
	findCenters(nr);
	
	// select from the candidates the one (if any) which is
	// enclosing and supporting:
	if (!selectBall(nr))
	  return false;
      }
      
      return true;
    }
    
  } // namespace Min_sphere_of_spheres_impl

} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_SUPPORT_SET_C
