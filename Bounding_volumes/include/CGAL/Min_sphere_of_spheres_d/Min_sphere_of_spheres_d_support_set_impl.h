// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Kaspar Fischer


#ifndef CGAL_MINIBALL_SUPPORTSET_C
#define CGAL_MINIBALL_SUPPORTSET_C

#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_support_set.h>

namespace CGAL_MINIBALL_NAMESPACE {

  namespace Min_sphere_of_spheres_d_impl {

    template<typename FT>
    inline bool reject(const FT& alpha,const FT& prev,
                       const Tag_false /* is_exact */) {
      using namespace Min_sphere_of_spheres_d_impl;
      return alpha < SqrOfEps<FT>::result() * sqr(prev);
    }

    template<typename FT,typename Pair>
    inline bool reject(const FT& alpha,const Pair&,
                       const Tag_true /* is_exact */) {
      return alpha == FT(0);
    }


    template<class Traits>
    bool Support_set<Traits>::push(const Sphere& ball) {
      CGAL_MINIBALL_DO_DEBUG(is_spanning_was_called = false);

      if ((std::ptrdiff_t)m > D)
        return false;

      b[m] = &ball;

      if (m == 0) {
        for (int j=0; j<D; ++j)
          d[0][j] = e[0][j] = f[0][j] = 0;
        sigma[0]     = 0;
        chi[0]       = -2;
        psi[0]       = FT(+4) * t.radius(ball);
        omega[0]     = FT(-2) * sqr(t.radius(ball));
        sol[1]       = t.radius(ball);
        discrim[1]   = 0;
        maxradius[0] = t.radius(ball);

      } else {
        // calculate $C_m$, storing it (temporarily) in u[m]:
        CIt c   = t.center_cartesian_begin(*b[m]),
            c_0 = t.center_cartesian_begin(*b[0]);
        for (int j=0; j<D; ++j) {
          u[m][j] = *c - *c_0;
          ++c; ++c_0;
        }

        // compute $\tau_{im}$ for $1<=i<m$:
        for (unsigned int i=1; i<m; ++i) {
          tau[i][m] = 0;
          for (int j=0; j<D; ++j)
            tau[i][m] += u[i][j]*u[m][j];
          tau[i][m] *= FT(2);
          tau[i][m] /= alpha[i];
        }

        // compute maxradius[m]:
        maxradius[m] = (std::max)(maxradius[m-1],t.radius(*b[m]));

        // calculate delta[m], eps[m] and phi[m] (by definition):
        const FT t1 = t.radius(*b[0]) - t.radius(*b[m]),
                     t2 = t.radius(*b[0]) + t.radius(*b[m]);
        phi[m]   = eps[m] = 0;
        delta[m] = -sigma[m-1];
        for (int j=0; j<D; ++j) {
          eps[m]   -= u[m][j]*e[m-1][j];
          phi[m]   -= u[m][j]*f[m-1][j];
          delta[m] += sqr(u[m][j]-d[m-1][j]);
        }
        phi[m] = FT(2)*(phi[m] - t1);
        eps[m] = t1*t2+FT(2)*eps[m];

        // fix u[m] to be $C_m-\q{C}_m$:
        // (This is only necessary for m>1 because $\q{C}_1=0$.)
        for (unsigned int i=1; i<m; ++i)
          for (int j=0; j<D; ++j)
            u[m][j] -= tau[i][m]*u[i][j];

        // calculate alpha[m]:
        alpha[m] = 0;
        for (int j=0; j<D; ++j)
          alpha[m] += sqr(u[m][j]);
        alpha[m] *= FT(2);

        // reject push if alpha[m] is to small:
        if (reject(alpha[m],sol[m],Is_exact()))
          return false;
        // calculate d[m], e[m] and f[m]:
        const FT da = delta[m]/alpha[m],
                     ea = eps[m]/alpha[m],
                     fa = phi[m]/alpha[m];
        for (int j=0; j<D; ++j) {
          d[m][j] = d[m-1][j] + da*u[m][j];
          e[m][j] = e[m-1][j] + ea*u[m][j];
          f[m][j] = f[m-1][j] + fa*u[m][j];
        }

        // compute sigma[m], chi[m], psi[m] and omega[m]:
        const FT de = delta[m]+eps[m];
        sigma[m] = sigma[m-1] + delta[m]*da/FT(2);
        chi[m]   = chi[m-1]   + phi[m]*fa;
        psi[m]   = psi[m-1]   + FT(2)*de*fa;
        omega[m] = omega[m-1] + de*de/alpha[m];

        // compute sol[m+1]:
        if (find_radius(Is_exact()) == false)
          return false;

      }

      ++m;
      return true;
    }

    template<class Traits>
    bool Support_set<Traits>::is_spanning() {
      CGAL_MINIBALL_DO_DEBUG(is_spanning_was_called = true);
      Result beta[D+1];

      // make sure at least one ball got pushed:
      CGAL_MINIBALL_ASSERT(m > 0);

      copy_n<D>(t.center_cartesian_begin(*b[0]),center);

      if (m > 1) {
        // compute the coeffients beta[i] and the center:
        for(unsigned int i=1; i<m; ++i) {
          beta[i] = (delta[i]+eps[i]+sol[m]*phi[i])/alpha[i];
          for (int j=0; j<D; ++j)
            center[j] += beta[i]*u[i][j];
        }

        // check whether the ball with center center and
        // radius sol[m] coincides with the miniball of the
        // pushed balls:
        Result gamma[D+1];
        Result mingamma(0);
        Result gamma0(1);

        for (unsigned int i=m-1; i>0; --i) {
          gamma[i] = beta[i];
          for (unsigned int j=i+1; j<m; ++j)
            gamma[i] -= gamma[j]*tau[i][j];
          gamma0 -= gamma[i];
          if (is_neg(gamma[i]-mingamma,discrim[m]))
            mingamma = gamma[i];
        }
        if (is_neg(gamma0-mingamma,discrim[m]))
          mingamma = gamma0;

        return !is_neg(mingamma,discrim[m]);
      }

      return true;
    }

  } // namespace Min_sphere_of_spheres_d_impl
} // namespace CGAL_MINIBALL_NAMESPACE

#endif // CGAL_MINIBALL_SUPPORTSET_CC
