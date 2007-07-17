// Copyright (c) 2007 Foundation for Research and Technology-Hellas (Greece).
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

CGAL_BEGIN_NAMESPACE

CGAL_APOLLONIUS_GRAPH_2_BEGIN_NAMESPACE

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
  typedef Arity_tag<3>             Arity;
  typedef Site_2                   argument_type;

public:
  inline
  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3) const
  {
    return Kernel().orientation_2_object()(s1.point(), s2.point(),
					   s3.point());
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

    FT Exp = det2x2_by_formula(xj, pj, xk, pk);
    FT Eyp = det2x2_by_formula(yj, pj, yk, pk);
    FT Erp = det2x2_by_formula(rj, pj, rk, pk);

    FT Exy = det2x2_by_formula(xj, yj, xk, yk);
    FT Exr = det2x2_by_formula(xj, rj, xk, rk);
    FT Eyr = det2x2_by_formula(yj, rj, yk, rk);

    FT Exy2 = 2 * det2x2_by_formula(xl, yl, xm, ym);

#if 0
    FT A = (Exp * Exr + Eyp * Eyr) * Exy2 + (Eyp * dx - Exp * dy) * Erp;
    FT B = Exy * Exy2 - Exp * dx - Eyp * dy;
    FT C = CGAL::square(Exp) + CGAL::square(Eyp) - CGAL::square(Erp);

    return sign_a_plus_b_x_sqrt_c(A, B, C);
#else
    Sign sA = CGAL::sign((Exp * Exr + Eyp * Eyr) * Exy2
			 + (Eyp * dx - Exp * dy) * Erp);

    FT B = Exy * Exy2 - Exp * dx - Eyp * dy;
    Sign sB = CGAL::sign(B);

    if ( sA == CGAL::ZERO ) { return sB; }
    if ( sB == CGAL::ZERO ) { return sA; }
    if ( sA == sB ) { return sA; }

    Sign s = CGAL::sign(CGAL::square(Exy2 * Exr - Erp * dy)
			+ CGAL::square(Exy2 * Eyr + Erp * dx)
			- CGAL::square(B));
    return sA * s;
#endif
  }

  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3, const Site_2& p1,
			 const Site_2& p2) const
  {
    Orientation o = predicate(s1, s2, s3, p1, p2);
    Orientation o_old = Base::operator()(s1, s2, s3, p1, p2);

    std::cerr << "Orientation predicate called" << std::endl;

    CGAL_assertion( o == o_old );
    return o;
  }

};

//--------------------------------------------------------------------

CGAL_APOLLONIUS_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_ORIENTATION8_C2_H
