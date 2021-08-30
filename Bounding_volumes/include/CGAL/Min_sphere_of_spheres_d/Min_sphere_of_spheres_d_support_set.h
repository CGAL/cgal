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

#ifndef CGAL_MINIBALL_SUPPORTSET
#define CGAL_MINIBALL_SUPPORTSET

#include <CGAL/license/Bounding_volumes.h>


#include <utility>
#include <functional>
#include <cmath>
#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_configure.h>

namespace CGAL_MINIBALL_NAMESPACE {
  namespace Min_sphere_of_spheres_d_impl {

    // The package provides specializations of (some of) the
    // algorithms for float and double, while the generic version
    // works with exact number types (which have no rounding
    // errors). So depending on the number type FT of the input (i.e.,
    // of the input point coordinates and radii) we choose a different
    // algorithm. In fact, not only the algorithm depends on FT but
    // also the type of the result (i.e., of the radius and the center
    // coordinates of the miniball): if FT is an exact number type,
    // the result is a Pair<FT> (because the radius and the center
    // coordinates may be degree-2 numbers which we encode as pairs);
    // if on the other hand FT is float or double, the result type
    // coincides with FT.
    //
    // The following struct Selector<FT> contains the appropriate
    // result type, and a flag telling whether the exact or floating-
    // point specialization of the chosen algorithm is to be run.
    template<typename FT>
    struct Selector {
      typedef Pair<FT> Result;
      typedef Tag_true Is_exact;
    };

    template<>
      struct Selector<float> {
        typedef float Result;
        typedef Tag_false Is_exact;
      };

    template<>
    struct Selector<double> {
      typedef double Result;
      typedef Tag_false Is_exact;
    };

    template<typename FT>
    struct Subtract_and_square {
      inline FT operator()(const FT x,const FT y) const {
        return sqr(x-y);
      }
    };

    template<typename FT>
    class Subtract_and_square_pair {
    private:
      const FT& d;

    public:
      Subtract_and_square_pair(const FT& d) : d(d) {}

      inline Pair<FT> operator()(const Pair<FT>& x,
                                     const FT& y) const {
        const FT t = x.first - y;
        return Pair<FT>(sqr(t)+sqr(x.second)*d,FT(2)*t*x.second);
      }
    };

    template<typename FT>
    struct Pair_to_double {
      const double root;

      Pair_to_double(const FT& disc) :
        root(std::sqrt(CGAL_MINIBALL_NTS to_double(disc))) {}

      inline double operator()(const Pair<FT>& p) const {
        return CGAL_MINIBALL_NTS to_double(p.first) +
               CGAL_MINIBALL_NTS to_double(p.second) * root;
      }
    };

    template<typename FT>
    struct Subtract_and_square_to_double {
      inline double operator()(const double x,const FT& y) const {
        return sqr(x - CGAL_MINIBALL_NTS to_double(y));
      }
    };

    template<class Traits>
    class Support_set {
    private: // some short hands:
      typedef typename Traits::FT FT;
      typedef typename Min_sphere_of_spheres_d_impl::
                       Selector<FT>::Result Result;
      typedef typename Min_sphere_of_spheres_d_impl::
                       Selector<FT>::Is_exact Is_exact;
      typedef typename Traits::Use_square_roots Use_sqrt;
      typedef typename Traits::Sphere Sphere;
      static const int D = Traits::D;
      typedef typename Traits::Cartesian_const_iterator CIt;

    public: // constructor:
      inline Support_set(Traits& traits);

    public: // access:
      inline const Result& radius() const;
      inline const Result *begin() const;
      inline const FT& disc() const;

    public: // containment test:
      template<typename InputIterator>
      bool contains(InputIterator c,const FT& r,
                    const FT tol,
                    const Tag_false /* is_exact */) const {
        // scale ball:
        const FT r1 = sol[m]*tol;

        // check radii:
        if (r > r1)
          return false;

        // compute the (squared) distance from center to c:
        const FT dist = inner_product_n<D>(center,c,
          FT(0),std::plus<FT>(),Subtract_and_square<FT>());

        // check containment:
        return dist <= sqr(r1-r);
      }

      template<typename InputIterator>
      bool contains(InputIterator c,const FT& r,
                    const double,const Tag_true /* is_exact */) const {
        typedef Pair<FT> P;

        // check radii:
        const P rd = sol[m]-r;
        if (is_neg(rd,discrim[m]))
          return false;

        // compute the (squared) distance from center to c:
        const P dist = inner_product_n<D>(center,c,
          P(0,0),std::plus<P>(),Subtract_and_square_pair<FT>(discrim[m]));

        // compute the square of rd:
        const P sqrRd(sqr(rd.first)+sqr(rd.second)*discrim[m],
                      FT(2)*rd.first*rd.second);

        // check containment:
        return is_neg_or_zero(dist-sqrRd,discrim[m]);
      }


    public: // modification:
      void reset();
      bool pivot(std::vector<const typename Traits::Sphere *>& l,
                 int& e,const int d);

    private: // modification:
      bool push(const Sphere& ball);
      inline void pop();
      bool is_spanning();

    private: // utility:
      inline bool find_radius(const Tag_false is_exact);
      inline bool find_radius(const Tag_true is_exact);

    private: // traits class:
      Traits& t;

    private: // for internal consisteny checks:
      #ifdef CGAL_MINIBALL_DEBUG
      // The following variable is true if and only if no ball has been
      // pushed so far, or is_spanning() has been called at least once and
      // neither pop() nor push() have been called since last is_spanning().
      bool is_spanning_was_called;
      #endif

    private: // member fields:
      unsigned int m;           // number of pushed balls
      const Sphere *b[D+1];     // pointers to pushed balls
      Result center[D+1];       // contains, when is_spanning() returns true,
                                // the center of the miniball

      // variables of the device:
      FT u[D+1][D];
      FT d[D+1][D];
      FT e[D+1][D];
      FT f[D+1][D];
      FT alpha[D+1];
      FT delta[D+1];
      FT eps[D+1];
      FT phi[D+1];
      FT sigma[D+1];
      FT chi[D+1];
      FT psi[D+1];
      FT omega[D+1];
      FT tau[D][D+1];
      Result sol[D+2];
      FT discrim[D+2];
      FT maxradius[D+1];
    };


    template<class Traits>
    Support_set<Traits>::Support_set(Traits& traits) : t(traits) {
      reset();
    }

    template<class Traits>
    void Support_set<Traits>::reset() {
      m = 0;
      sol[0] = FT(-1);
      discrim[0] = FT(0);
      CGAL_MINIBALL_DO_DEBUG(is_spanning_was_called = true);
    }

    template<class Traits>
    const typename Support_set<Traits>::Result&
      Support_set<Traits>::radius() const {
      CGAL_MINIBALL_ASSERT(is_spanning_was_called);
      return sol[m];
    }

    template<class Traits>
    const typename Support_set<Traits>::Result
      *Support_set<Traits>::begin() const {
      CGAL_MINIBALL_ASSERT(is_spanning_was_called);
      return center;
    }

    template<class Traits>
    const typename Support_set<Traits>::FT&
      Support_set<Traits>::disc() const {
      CGAL_MINIBALL_ASSERT(is_spanning_was_called);
      return discrim[m];
    }

    template<class Traits>
    void Support_set<Traits>::pop() {
      CGAL_MINIBALL_ASSERT(m>0);
      CGAL_MINIBALL_DO_DEBUG(is_spanning_was_called = false);
      --m;
    }

    template<class Traits>
    bool Support_set<Traits>::find_radius(const Tag_false /* is_exact */) {
      using namespace Min_sphere_of_spheres_d_impl;

      // find discriminant:
      discrim[m+1] = sqr(psi[m]) - 4*chi[m]*omega[m];
      if (discrim[m+1] < Min_float)
        return false;

      // handle degenerate quadratic:
      if (std::abs(chi[m]) < Min_float) {
        if (std::abs(psi[m]) < Min_float)
          return false;
        sol[m+1] = -omega[m]/psi[m];
        return sol[m+1]*Tol<FT>::result() >= maxradius[m];
      }

      // calculate the two solutions:
      FT sd = std::sqrt(discrim[m+1]);
      if (psi[m] > 0)
        sd = -sd;
      Result sols[2] = { (sd-psi[m])     / (FT(2)*chi[m]),
                          FT(2)*omega[m] / (sd-psi[m]) };

      // order the solutions (*):
      if (sols[1] < sols[0])
        std::swap(sols[0],sols[1]);

      // eliminate negative solutions:
      if (sols[0]*Tol<FT>::result() >= maxradius[m]) {
        sol[m+1] = sols[0];
        return true;
      }
      sol[m+1] = sols[1];
      return sol[m+1]*Tol<FT>::result() >= maxradius[m];
    }

    template<class Traits>
    bool Support_set<Traits>::find_radius(const Tag_true /* is_exact */) {
      // find discriminant:
      discrim[m+1] = sqr(psi[m]) - FT(4)*chi[m]*omega[m];
      if (discrim[m+1] < FT(0))
        return false;

      // handle degenerate quadratic:
      if (chi[m] == FT(0)) {
        if (psi[m] == FT(0))
          return false;
        sol[m+1] = Pair<FT>(-omega[m]/psi[m],FT(0));
        return sol[m+1].first >= maxradius[m];
      }

      // calculate the two solutions:
      const FT tmp = FT(-1)/(FT(2)*chi[m]);
      Result sols[2] = {
        Pair<FT>(psi[m]*tmp, tmp),
        Pair<FT>(psi[m]*tmp,-tmp)
      };

      // order the solutions (*):
      if (chi[m] < FT(0))
        std::swap(sols[0],sols[1]);
      CGAL_MINIBALL_ASSERT(is_neg_or_zero(sols[0]-sols[1],discrim[m+1]));

      // eliminate negative solutions:
      if (!is_neg(sols[0]-maxradius[m],discrim[m+1])) {
        sol[m+1] = sols[0];
        return true;
      }
      sol[m+1] = sols[1];
      return !is_neg(sols[1]-maxradius[m],discrim[m+1]);
    }

  } // namespace Min_sphere_of_spheres_d_impl
} // namespace CGAL_MINIBALL_NAMESPACE

#include <CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_support_set_impl.h>

#endif // CGAL_MINIBALL_SUPPORTSET
