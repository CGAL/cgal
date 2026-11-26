// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

#ifndef CGAL_HYPERBOLIC_SURFACE_DELAUNAY_TRAITS_2
#define CGAL_HYPERBOLIC_SURFACE_DELAUNAY_TRAITS_2

#ifdef CGAL_USE_GMP

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Complex_number.h>
#include <CGAL/Root_of_traits.h>
#include <boost/numeric/interval.hpp>
#include <CGAL/Gmpfr.h>
//#include <CGAL/Gmpq.h>

namespace CGAL {

  /* returns an approximation of the cc of the triangle pqr, which coordinates
 are approximations at relative precision \f$ 2^{-53\times p}\f$ of the exact
 circumcenter.  Requirements: FT==Gmpq
  */
  template <typename Traits>
    class  Construct_approximate_hyperbolic_circumcenter_2_f {
  public:
    typedef typename Traits::FT                          FT;
    typedef typename Traits::Hyperbolic_point_2          Hyperbolic_point_2;
    typedef typename Traits::Hyperbolic_Voronoi_point_2          Hyperbolic_Voronoi_point_2;
    typedef Complex_number<FT>                                          Complex;

  Construct_approximate_hyperbolic_circumcenter_2_f(Traits gt = Traits()): _gt(gt) {}
    Hyperbolic_point_2 operator()(const Hyperbolic_point_2& p,
                                  const Hyperbolic_point_2& q,
                                  const Hyperbolic_point_2& r) const
    {
      typename Traits::Construct_hyperbolic_circumcenter_2 chc = _gt.construct_hyperbolic_circumcenter_2_object();
      auto c = chc(p,q,r);


  if constexpr(std::is_same_v<FT, Gmpq>) {
      FT x;
      FT y;
      unsigned p;
      p = _gt.approximation_precision() * _gt.DOUBLE_PREC;
      Gmpfr a_0 = Gmpfr( c.x().a0().numerator(), p) / Gmpfr( c.x().a0().denominator(), p);
      Gmpfr a_1 = Gmpfr( c.x().a1().numerator(), p) / Gmpfr( c.x().a1().denominator(), p);
      Gmpfr r = Gmpfr( c.x().root().numerator(), p) / Gmpfr( c.x().root().denominator(), p);
      x = a_0 + a_1 * sqrt(r);
      a_0 = Gmpfr( c.y().a0().numerator(), p) / Gmpfr( c.y().a0().denominator(), p);
      a_1 = Gmpfr( c.y().a1().numerator(), p) / Gmpfr( c.y().a1().denominator(), p);
      r = Gmpfr( c.y().root().numerator(), p) / Gmpfr( c.y().root().denominator(), p);
      y = a_0 + a_1 * sqrt(r);
      if (!(norm(Complex(x, y)) < FT(1))) {std::cout<< "WARNING: THE CONSTRUCTION OF THE CIRCUMCENTER FAILED, THE APPROXIMATION PRECISION SHOULD BE INCREASED\n";}
      Hyperbolic_point_2 c_approx =  Hyperbolic_point_2(x, y);
      return c_approx;
    }
  else
    {
      Hyperbolic_point_2 c_approx;
      return c_approx(FT(to_double(c.x())), FT(to_double(c.y())));
    }
    }
  private:
    Traits _gt;
  }; // end class  Construct_approximate_hyperbolic_circumcenter_2_fct


template<class HyperbolicSurfaceTraitsClass>
class Hyperbolic_surface_Delaunay_traits_2
: public HyperbolicSurfaceTraitsClass
{
 public:
  typedef typename HyperbolicSurfaceTraitsClass::FT                          FT;
  typedef typename HyperbolicSurfaceTraitsClass::Hyperbolic_point_2          Hyperbolic_point_2;
  typedef typename HyperbolicSurfaceTraitsClass::Hyperbolic_Voronoi_point_2          Hyperbolic_Voronoi_point_2;
  typedef Complex_number<FT>                                          Complex;
  typedef Hyperbolic_surface_Delaunay_traits_2<HyperbolicSurfaceTraitsClass>    Self;
  typedef Construct_approximate_hyperbolic_circumcenter_2_f<Self>   Construct_approximate_hyperbolic_circumcenter_2;

  static const unsigned DOUBLE_PREC = 53;

  Construct_approximate_hyperbolic_circumcenter_2
    construct_approximate_hyperbolic_circumcenter_2_object() const
  { return Construct_approximate_hyperbolic_circumcenter_2(*this); };

 Hyperbolic_surface_Delaunay_traits_2() : approx_precision(1) {};

  unsigned approximation_precision() const {return approx_precision;}
  void approximation_precision(unsigned p) {approx_precision = p;}

 private:
  unsigned approx_precision;

 };

} // namespace CGAL

#endif //CGAL_USE_GMP

#endif // CGAL_HYPERBOLIC_SURFACE_DELAUNAY_TRAITS_2
