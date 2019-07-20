// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_HYPERBOLA_SEGMENT_2_H
#define CGAL_HYPERBOLA_SEGMENT_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Hyperbola_2.h>

namespace CGAL {

template < class Gt >
class Hyperbola_segment_2 : public Hyperbola_2< Gt >
{
public:
  typedef CGAL::Hyperbola_2<Gt>               Base;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2                Segment_2;
  typedef typename Base::FT                       FT;

  using Base::t;
  using Base::f;

#if 0
  typedef CGAL::Hyperbola_2<Point,Weight>         Base;
  typedef typename Base::Weighted_point           Weighted_point;
#endif
  //  typedef typename R::RT          FT;
  //  typedef double                                  FT;
  //  typedef CGAL::Point_2< Cartesian<double> >      Point_2;
  //  typedef CGAL::Segment_2< Cartesian<double> >    Segment_2;

protected:
  Point_2 p1, p2;

  template< class Stream >
  inline
  void draw_line(Stream &W) const
  {
#if 0
    FT s[2];

    s[0] = t(p1);
    s[1] = t(p2);

    Point_2 p[2];
    for (int i = 0; i < 2; i++)  p[i] = f(s[i]);

    W << Segment_2(p[0], p[1]);
#else
    W << Segment_2(p1, p2);
#endif
  }

  inline
  Point_2 midpoint() const
  {
    return Hyperbola_2< Gt >::midpoint(p1, p2);
  }

public:
  Hyperbola_segment_2() : Hyperbola_2< Gt >() {}

  Hyperbola_segment_2(const Site_2 &f1,	const Site_2 &f2,
		      const Point_2 &p1, const Point_2 &p2)
    : Hyperbola_2< Gt >(f1, f2)
  {
    this->p1 = p1;
    this->p2 = p2;
  }

  void generate_points(std::vector<Point_2>& p) const
  {
    if ( CGAL::is_zero(this->r) ) {
      p.push_back(p1);
      p.push_back(p2);
      return;
    }


    //    FT STEP = W.width() / 100.0;

    FT s0, s1;

    s0 = t(p1);
    s1 = t(p2);

    if (CGAL::compare(s0, s1) == LARGER) {
      std::swap(s0, s1);
    }

    p.clear();

    if ( !(CGAL::is_positive(s0)) &&
	 !(CGAL::is_negative(s1)) ) {
      FT tt;
      int k;

      p.push_back( this->o );
      k = 1;
      tt = FT(-this->STEP);
      while ( CGAL::compare(tt, s0) == LARGER ) {
	p.insert( p.begin(), f(tt) );
	k--;
	tt = -FT(k * k) * this->STEP;
      }
      p.insert( p.begin(), f(s0) );

      k = 1;
      tt = FT(this->STEP);
      while ( CGAL::compare(tt, s1) == SMALLER ) {
	p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * this->STEP;
      }
      p.push_back( f(s1) );
    } else if ( !(CGAL::is_negative(s0)) &&
		!(CGAL::is_negative(s1)) ) {
      FT tt;
      int k;


      p.push_back( f(s0) );

      tt = s0;
      k = int(CGAL::to_double(CGAL::sqrt(tt / this->STEP)));

      while ( CGAL::compare(tt, s1) == SMALLER ) {
	if ( CGAL::compare(tt, s0) != SMALLER )
	  p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * this->STEP;
      }
      p.push_back( f(s1) );
    } else {
      FT tt;
      int k;

      p.push_back( f(s1) );

      tt = s1;
      k = int(CGAL::to_double(-CGAL::sqrt(-tt / this->STEP)));

      while ( CGAL::compare(tt, s0) == LARGER ) {
	if ( CGAL::compare(tt, s1) != LARGER )
	  p.push_back( f(tt) );
	k--;
	tt = -FT(k * k) * this->STEP;
      }
      p.push_back( f(s0) );
    }
  }

  template< class Stream >
  void draw(Stream &W) const
  {
    std::vector<Point_2> p;
    generate_points(p);

    for (unsigned int i = 0; i < p.size() - 1; i++) {
      W << Segment_2(p[i], p[i+1]);
    }
  }


};


template< class Stream, class Gt >
inline
Stream& operator<<(Stream &s, const Hyperbola_segment_2<Gt>& H)
{
  H.draw(s);
  return s;
}

} //namespace CGAL

#endif // CGAL_HYPERBOLA_SEGMENT_2_H
