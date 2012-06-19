// Copyright (c) 2007 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_GRAPH_2_ORIENTATION8_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_ORIENTATION8_C2_H

#include <CGAL/determinant.h>
#include <CGAL/Apollonius_graph_2/Orientation_2.h>

//--------------------------------------------------------------------

namespace CGAL {

namespace ApolloniusGraph_2 {

template<class K, class MTag>
class Orientation8_C2
  : public Orientation_2<K,MTag>
{
private:
  typedef Orientation_2<K,MTag>    Base;

public:
  typedef K                        Kernel;
  typedef MTag                     Method_tag;
  typedef typename K::Site_2       Site_2;
  typedef typename K::Point_2      Point_2;
  typedef typename K::Orientation  Orientation;
  typedef typename K::FT           FT;

  typedef Orientation              result_type;
  typedef Site_2                   argument_type;

public:
  inline
  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3) const
  {
    return Kernel().orientation_2_object()(s1.point(), s2.point(),
					   s3.point());
  }

  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B,
		     const FT& Exp, const FT& Eyp, const FT& Erp,
		     const FT& Exy2, const FT& Exr, const FT& Eyr,
		     const FT dx, const FT& dy,
		     const Field_with_sqrt_tag&) const
  {
    FT G = CGAL::square(Exp) + CGAL::square(Eyp) - CGAL::square(Erp);
    return CGAL::sign(A + B * CGAL::sqrt(G));
  }

  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B,
		     const FT& Exp, const FT& Eyp, const FT& Erp,
		     const FT& Exy2, const FT& Exr, const FT& Eyr,
		     const FT dx, const FT& dy,
		     const Integral_domain_without_division_tag&) const
  {
    Sign sA = CGAL::sign(A);
    Sign sB = CGAL::sign(B);

    if ( sA == CGAL::ZERO ) { return sB; }
    if ( sB == CGAL::ZERO ) { return sA; }
    if ( sA == sB ) { return sA; }

    Sign s = CGAL::sign(CGAL::square(Exy2 * Exr - Erp * dy)
			+ CGAL::square(Exy2 * Eyr + Erp * dx)
			- CGAL::square(B));
    return sA * s;
  }

  Orientation predicate(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3, const Site_2& p1,
			 const Site_2& p2) const
  {
    // computes the orientation of the Voronoi vertex of s1, s2, s3 and
    // the points p1 and p2
    FT xj = s2.x() - s1.x();
    FT xk = s3.x() - s1.x();

    FT xl = p1.x() - s1.x();
    FT xm = p2.x() - s1.x();


    FT yj = s2.y() - s1.y();
    FT yk = s3.y() - s1.y();

    FT yl = p1.y() - s1.y();
    FT ym = p2.y() - s1.y();


    FT dx = xl - xm;
    FT dy = yl - ym;


    FT rj = s2.weight() - s1.weight();
    FT rk = s3.weight() - s1.weight();

    FT pj = CGAL::square(xj) + CGAL::square(yj) - CGAL::square(rj);
    FT pk = CGAL::square(xk) + CGAL::square(yk) - CGAL::square(rk);

    FT Exp = determinant(xj, pj, xk, pk);
    FT Eyp = determinant(yj, pj, yk, pk);
    FT Erp = determinant(rj, pj, rk, pk);

    FT Exy = determinant(xj, yj, xk, yk);
    FT Exr = determinant(xj, rj, xk, rk);
    FT Eyr = determinant(yj, rj, yk, rk);

    FT Exy2 = 2 * determinant(xl, yl, xm, ym);

    FT A = (Exp * Exr + Eyp * Eyr) * Exy2 + (Eyp * dx - Exp * dy) * Erp;
    FT B = Exy * Exy2 - Exp * dx - Eyp * dy;

    return sqrt_ext_sign(A, B, Exp, Eyp, Erp,
			 Exy2, Exr, Eyr, dx, dy, Method_tag());
  }

  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3, const Site_2& p1,
			 const Site_2& p2) const
  {
    Orientation o = predicate(s1, s2, s3, p1, p2);
#ifndef NDEBUG
    Orientation o_old = Base::operator()(s1, s2, s3, p1, p2);

    CGAL_assertion( o == o_old );
#endif
    return o;
  }

};

//--------------------------------------------------------------------

template<class K, class MTag>
class Constructive_orientation8_C2
{
public:
  typedef K                        Kernel;
  typedef MTag                     Method_tag;
  typedef typename K::Site_2       Site_2;
  typedef typename K::Point_2      Point_2;
  typedef typename K::Orientation  Orientation;
  typedef typename K::FT           FT;

  typedef Orientation              result_type;
  typedef Site_2                   argument_type;

private:
  FT s1x, s1y;
  FT xj, xk, yj, yk, rj, rk, nj, nk, pj, pk, Exp, Eyp, Erp, Exy, Exr, Eyr, A1;
  Orientation o_sym;

public:
  Constructive_orientation8_C2(const Site_2& s1, const Site_2& s2,
			       const Site_2& s3, bool use_xj)
  {
    s1x = s1.x();
    s1y = s1.y();
    
    xj = s2.x() - s1.x();
    xk = s3.x() - s1.x();

    yj = s2.y() - s1.y();
    yk = s3.y() - s1.y();

    rj = s2.weight() - s1.weight();
    rk = s3.weight() - s1.weight();

    nj = CGAL::square(xj) + CGAL::square(yj);
    nk = CGAL::square(xk) + CGAL::square(yk);

    pj = nj - CGAL::square(rj);
    pk = nk - CGAL::square(rk);

    Exp = determinant(xj, pj, xk, pk);
    Eyp = determinant(yj, pj, yk, pk);
    Erp = determinant(rj, pj, rk, pk);

    Exy = determinant(xj, yj, xk, yk);
    Exr = determinant(xj, rj, xk, rk);
    Eyr = determinant(yj, rj, yk, rk);

    A1 = Exp * Exr + Eyp * Eyr;

    FT A, B, norm;
    if ( use_xj ) {
      A = (-Eyp * xj + Exp * yj) * Erp;
      B = Exp * xj + Eyp * yj;
      norm = nj;
    } else {
      A = (-Eyp * xk + Exp * yk) * Erp;
      B = Exp * xk + Eyp * yk;
      norm = nk;
    }

    o_sym = sqrt_ext_sign(A, B, norm, Method_tag());
  }

private:
  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B, const FT& Exy2, 
		     const FT& dx, const FT& dy,
		     const Field_with_sqrt_tag&) const
  {
    FT G = CGAL::square(Exp) + CGAL::square(Eyp) - CGAL::square(Erp);
    return CGAL::sign(A + B * CGAL::sqrt(G));
  }

  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B, const FT&,
		     const Field_with_sqrt_tag&) const
  {
    FT G = CGAL::square(Exp) + CGAL::square(Eyp) - CGAL::square(Erp);
    return CGAL::sign(A + B * CGAL::sqrt(G));
  }


  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B, const FT& norm,
		     const Integral_domain_without_division_tag&) const
  {
    Sign sA = CGAL::sign(A);
    Sign sB = CGAL::sign(B);

    if ( sA == CGAL::ZERO ) { return sB; }
    if ( sB == CGAL::ZERO ) { return sA; }
    if ( sA == sB ) { return sA; }

    Sign s = CGAL::sign(CGAL::square(Erp) * norm - CGAL::square(B));
    return sA * s;
  }


  inline
  Sign sqrt_ext_sign(const FT& A, const FT& B, const FT& Exy2,
		     const FT dx, const FT& dy,
		     const Integral_domain_without_division_tag&) const
  {
    Sign sA = CGAL::sign(A);
    Sign sB = CGAL::sign(B);

    if ( sA == CGAL::ZERO ) { return sB; }
    if ( sB == CGAL::ZERO ) { return sA; }
    if ( sA == sB ) { return sA; }

    Sign s = CGAL::sign(CGAL::square(Exy2 * Exr - Erp * dy)
			+ CGAL::square(Exy2 * Eyr + Erp * dx)
			- CGAL::square(B));
    return sA * s;
  }

  Orientation predicate(const Site_2& p1, const Site_2& p2) const
  {
    // computes the orientation of the Voronoi vertex of s1, s2, s3 and
    // the points p1 and p2
    FT xl = p1.x() - s1x;
    FT xm = p2.x() - s1x;

    FT yl = p1.y() - s1y;
    FT ym = p2.y() - s1y;

    FT dx = xl - xm;
    FT dy = yl - ym;

    FT Exy2 = 2 * determinant(xl, yl, xm, ym);

    FT A = A1 * Exy2 + (Eyp * dx - Exp * dy) * Erp;
    FT B = Exy * Exy2 - Exp * dx - Eyp * dy;

    return sqrt_ext_sign(A, B, Exy2, dx, dy, Method_tag());
  }

public:
  inline
  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3) const
  {
    return Kernel().orientation_2_object()(s1.point(), s2.point(),
					   s3.point());
  }

  inline
  Orientation operator()(const Site_2& p1, const Site_2& p2) const
  {
    Orientation o = predicate(p1, p2);
    return o;
  }

  inline
  Orientation operator()() const {
    return o_sym;
  }
};

//--------------------------------------------------------------------


} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_ORIENTATION8_C2_H
