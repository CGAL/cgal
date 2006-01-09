// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_APPROXIMATE_MIN_ELLIPSOID_D_H
#define CGAL_APPROXIMATE_MIN_ELLIPSOID_D_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include <CGAL/Approximate_min_ellipsoid_d/Khachiyan_approximation.h>

namespace CGAL {

  template<class Traits>
  class Approximate_min_ellipsoid_d
  // An instance of class Approximate_min_ellipsoid_d represents an
  // eps-approximation of the smallest enclosing ellipsoid of the
  // point set P in R^d, that is, an ellipsoid E satisfying
  //
  //  (i) E contains P,
  //
  //  (ii) the volume of E is at most (1+eps) times larger than the
  //    volume of the smallest enclosing ellipsoid of P.
  {
  public: // types:
    typedef typename Traits::FT                       FT;
    typedef typename Traits::Point                    Point;
    typedef typename Traits::Cartesian_const_iterator Cartesian_const_iterator;

  private:
    typedef std::vector<Point>                        Point_list;
    typedef typename Traits::Cartesian_const_iterator C_it;
    
  protected: // member variables:
    Traits tco;                           // traits class object
    const Point_list P;                   // the input points in R^d
    int d;                                // dimension of the input points
    const double eps;                     // the desired epsilon
    double e_eps;                         // (see below)
    
    // We obtain our eps-approximating ellipsoid by embedding the input
    // points P into R^{d+1} by mapping p to (p,1).  Then we compute in
    // R^{d+1} an e_eps-approximating centrally symmetric ellipsoid E for
    // the embedded points from which the desired eps-approximating
    // ellipsoid of the original points can be obtained by projection.
    // The following variable E represents the e_eps-approximating
    // centrally symmetric ellipsoid in R^{d+1}:
    Khachiyan_approximation<true,Traits> *E;

    // When the input points do not affinely span the whole space
    // (i.e., if dim(aff(P)) < d), then the smallest enclosing
    // ellipsoid of P has no volume in R^d and so the points are
    // called "degnerate" (see is_degenerate()) below.

    // As discussed below (before (*)), the centrally symmetric ellipsoid
    // E':= sqrt{(1+a_eps)(d+1)} E contains (under exact arithmetic) the
    // embedded points. (Here, a_eps is the value obtained by
    // achieved_epsilon().) Denoting by M the defining matrix of E; we
    // then have
    //
    //      p^T M p <= 1
    // 
    // for all p in E.  Since this is equivalent to
    //
    //     sqrt{alpha}p^T alpha M sqrt{alpha}p <= 1,
    //
    // for alpha = 1/((1+a_eps)(d+1)), we see that E' has alpha M as its
    // defining matrix. Consequently, the routines defining_matrix(),
    // defining_vector(), and defining_scalar() will return numbers that
    // need to be scaled (by the user) with the factor alpha. (We do not
    // perform the scaling ourselves because we cannot do it exactly.)

  public: // construction & destruction:

    template<typename InputIterator>
    Approximate_min_ellipsoid_d(double eps,
			    InputIterator first,InputIterator last,
			    const Traits& traits = Traits())
    // Given a range [first,last) of n points P, constructs an
    // (1+eps)-approximation of the smallest enclosing ellipsoid of P.
      : tco(traits), P(first,last), eps(eps)
    {
      CGAL_APPEL_LOG("appel",
		     "Entering Approximate_min_ellpsoid_d." << std::endl);

      // fetch ambient dimension:
      d = tco.dimension(P[0]);

      // The ellipsoid E produced by Khachiyan's algorithm has the
      // property that E':= sqrt{(1+e_eps) (d+1)} E contains all
      // embedded points e_P and has volume bounded by
      //
      //   vol(E') <= (1+e_eps)^{(d+1)/2} vol(E_emb*),            (*)
      //
      // where E_emb* is the smallest centrally symmetric enclosing
      // ellipsoid of the embedded points e_P.  By requiring that (1+eps)
      // <= (1+e_eps)^{(d+1)/2}, we get
      //
      //   e_eps <= (1+eps)^{2/(d+1)} - 1                         (**)
      //
      // and with this e_eps, we have vol(E') <= (1+eps) vol(E_emb*).
      // An argument by Khachiyan (see "Rounding of polytopes in the
      // real number model of computation", eq. (4.1)) then guarantees
      // that the intersection E_int of E' with the hyperplane { (x,y)
      // in R^{d+1} | y = 1} is an ellipsoid satisfying vol(E_int) <=
      // (1+eps) vol(E*), where E* is the smallest enclosing ellipsoid
      // of the original points P.  According to (**) we thus set e_eps
      // to (a lower bound on) (1+eps)^{2/(d+1)} - 1:
      FPU_CW_t old = FPU_get_and_set_cw(CGAL_FE_TOWARDZERO); // round to zero
      e_eps = std::exp(2.0/(d+1)*std::log(1.0+eps))-1.0;
      FPU_set_cw(old);                                       // restore
							     // rounding mode
      
      // Find e (1+e_eps)-approximation for the embedded points.  This
      // only works when the points affinely span R^{d+1}.
      E = new Khachiyan_approximation<true,Traits>(d,P.size(),tco);
      const bool is_deg = !E->add(P.begin(),P.end(),e_eps);

      // debugging:
      CGAL_APPEL_ASSERT(is_deg == E->is_degenerate());
      CGAL_APPEL_LOG("appel",
		     "  Input points are " << (is_deg? "" : "not ") <<
		     "degnerate." << std::endl);

      if (is_deg)
	find_lower_dimensional_approximation();

      CGAL_APPEL_LOG("appel",
		     "Leaving Approximate_min_ellipsoid_d." << std::endl);
    }
    
    ~Approximate_min_ellipsoid_d()
    {
      if (E != static_cast<Khachiyan_approximation<true,Traits> *>(0))
	delete E;
    }

  public: // access:

    unsigned int number_of_points() const
    // Returns the number of points, i.e., |P|.
    {
      return P.size();
    }
    
    bool is_empty() const
    // Returns true iff the approximate ellipsoid is empty (which
    // implies degeneracy).  This is the case iff the number of input
    // points was zero at construction time.
    {
      return P.size() == 0;
    }
    
  private: // access:
    bool is_degenerate() const
    // Returns true iff the approximate ellipsoid is degenerate, i.e.,
    // iff the dimension of the affine hull of S doesn't match the
    // dimension of the ambient space.
    {
      return E->is_degenerate();
    }

  public: // access:

    bool is_full_dimensional() const
    // Returns !is_degenerate().
    {
      return !is_degenerate();
    }

    FT defining_matrix(int i,int j) const
    // Returns the entry M(i,j) of the symmetric matrix M in the
    // representation
    //
    //    E = { x | x^T M x + x^T m + mu <= 0 }
    //
    // of the computed approximation. More precisely, the routine does not
    // return M(i,j) but the number (1+achieved_epsilon())*(d+1)*M(i,j).
    //
    // Precondition: !is_degenerate() && 0<=i<d && 0<=j<d
    {
      
      CGAL_APPEL_ASSERT(!is_degenerate() && 0<=i && i<d && 0<=j && j<d);
      return E->matrix(i,j);
    }

    FT defining_vector(int i) const
    // Returns the entry m(i) of the vector m in the representation
    //
    //    E = { x | x^T M x + x^T m + mu <= 0 }
    //
    // of the computed approximation. More precisely, the routine does not
    // return m(i) but the number (1+achieved_epsilon())*(d+1)*m(i).
    //
    // Precondition: !is_degenerate() && 0<=i<d
    {
      CGAL_APPEL_ASSERT(!is_degenerate() && 0<=i && i<d);
      return FT(2)*E->matrix(d,i); // Note: if FT is double, the
				   // multiplication by 2.0 is exact.
    }
    
    FT defining_scalar() const
    // Returns the number mu in the representation
    //
    //    E = { x | x^T M x + x^T m + mu <= 0}
    //
    // of the computed approximation. More precisely, the routine does not
    // return mu but the number (1+achieved_epsilon())*(d+1)*mu+1.
    //
    // Precondition: !is_degenerate()
    {
      CGAL_APPEL_ASSERT(!is_degenerate());
      return E->matrix(d,d);
    }

    double achieved_epsilon() const
    // Returns the approximation ratio; more precisely, this returns a
    // number r such that the computed ellipsoid is a (1+r)-approximation
    // of MEL(P).
    //
    // Precondition: !is_degenerate()
    {
      CGAL_APPEL_ASSERT(!is_degenerate());
      
      // From (*) we known that Khachian's algorithm produces a
      // centrally-symmetric ellipsoid in R^{d+1} fulfilling
      //
      //  vol(E') <= (1+k_eps)^{(d+1)/2} vol(E_emb*),
      //
      // where k_eps = E->exact_epsilon(). The projecting argument
      // mentioned after (**) says that these approximation ratio also
      // applies to the projected (d-dimensional) ellipsoids, and so the
      // actual approximation ratio we obtain is
      //
      //   ratio:= (1+k_eps)^{(d+1)/2}.
      //
      // So all we need to do is compute the smallest (let's say: a small)
      // double number larger or equal to ratio.
      //
      // Todo: in the following implementation it would maybe be more
      // accurate to use CGAL::Interval_nt to compute sqrt(1+k_eps) and
      // then integer-exponentiate this with exponent d+1.
      const double k_eps = E->exact_epsilon();
      FPU_CW_t old = FPU_get_and_set_cw(CGAL_FE_UPWARD); // round up
      const double eps = std::exp((d+1)*0.5*std::log(1.0+k_eps))-1.0;
      FPU_set_cw(old);                                   // restore
      
      // Todo: following line should be uncommented as soon as
      // the bug with eps < 0 is fixed.
      // CGAL_APPEL_ASSERT(eps >= 0.0);
      return eps;
    }

    Traits traits() const
    {
      return tco;
    }

  public: // miscellaneous:
    
    bool is_valid(bool verbose) const
    // Returns true if and only if the computed ellipsoid is indeed an
    // approximate ellipsoid, that is ... Todo.
    {
      return E->is_valid(verbose);
    }
      
  public: // miscellaneous 3D support:

    typedef double *Axis_length_iterator_3;
    typedef double *Axis_direction_coordinates_iterator_3;
    
    Axis_length_iterator_3 axis_length_begin_3();
    // Returns an iterator of value type double to the first of the
    // three lengths of the three axes of the computed ellipsoid.
    //
    // Note: the lengths you obtain this way are approximations to the
    // real lengths and come without any guarantee on the relative
    // error.
    //
    // Precondition: d==3 && !is_degenerate().

    Axis_direction_coordinates_iterator_3 axis_direction_coordinates_begin_3();
    // Returns an iterator to the first of three vectors the coordilengths of the 3 main axes of the
    // computed ellipsoid.
    //
    // Precondition: d==3 && !is_degenerate().
    

  public: // internal debugging routines:

    void write_eps(const std::string& name) const;
    // Writes the and the point set P to an EPS file.  Returned is the
    // id under which the file was stored (filename 'id.eps').
    //
    // Precondition: d==2 && !is_degenerate().

  private: // internal routines:

    void find_lower_dimensional_approximation(); // (does nothing right now)
  };
  
}

#include <CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d.C>

#endif // CGAL_APPROXIMATE_MIN_ELLIPSOID_D_H
