// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_PARABOLA_SEGMENT_2_H
#define CGAL_PARABOLA_SEGMENT_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Parabola_2.h>

namespace CGAL {

namespace Qt {
  template <typename K> class PainterOstream;
}

template < class Gt >
class Parabola_segment_2 : public Parabola_2< Gt >
{
  typedef CGAL::Parabola_2<Gt>            Base;
  typedef typename Base::Site_2           Site_2;
  typedef typename Base::FT               FT;
  typedef typename Base::Point_2          Point_2;
  typedef typename Base::Segment_2        Segment_2;
  typedef typename Base::Line_2           Line_2;

  using Base::t;
  using Base::f;

protected:
  Point_2 p1, p2;

public:
  Parabola_segment_2() : Parabola_2< Gt >() {}

  template<class ApolloniusSite>
  Parabola_segment_2(const ApolloniusSite &p, const Line_2 &l,
                     const Point_2 &p1, const Point_2 &p2)
    : Parabola_2< Gt >(p, l)
  {
    this->p1 = p1;
    this->p2 = p2;
  }

  Parabola_segment_2(const Point_2 &p, const Line_2 &l,
                     const Point_2 &p1, const Point_2 &p2)
    : Parabola_2< Gt >(p, l)
  {
    this->p1 = p1;
    this->p2 = p2;
  }

  int compute_k(const FT tt, const FT STEP) const {
    return int(CGAL::to_double(CGAL::sqrt(tt / STEP)));
  }

  // s0 and s1 define a desired drawing "range"
  void generate_points(std::vector<Point_2>& p,
                       const FT STEP,
                       FT s0, FT s1) const
  {
    CGAL_precondition(STEP > 0);

    p.clear();

    if (CGAL::compare(s0, s1) == LARGER)
      std::swap(s0, s1);

    // This is a parabola segment that exists between only between p1 and p2 so we gotta crop
    // the desired range to actually fit the parabola segment
    FT tp1 = t(p1), tp2 = t(p2);

    if (CGAL::compare(tp1, tp2) == LARGER)
      std::swap(tp1, tp2);

    if(tp2 < s0 || s1 < tp1) // nothing to draw because the ranges are completely disjoint
      return;

    s0 = (std::max)(s0, tp1);
    s1 = (std::min)(s1, tp2);

    if ( !(CGAL::is_positive(s0)) && !(CGAL::is_negative(s1)) )
    {
      FT tt = - STEP;
      int k = -1;

      p.push_back( this->o );

      // no need to check tt < s1 since we have tt < 0, s1 >= 0 and tt is moving towards -inf
      while ( CGAL::compare(tt, s0) == LARGER )
      {
        p.insert( p.begin(), f(tt) );
        --k;
        tt = - FT(k * k) * STEP;
      }
      p.insert( p.begin(), f(s0) );

      k = 1;
      tt = STEP;

      // no need to check tt > s0 since we have tt > 0, s0 <= 0 and tt is moving towards +inf
      while ( CGAL::compare(tt, s1) == SMALLER )
      {
        p.push_back( f(tt) );
        ++k;
        tt = FT(k * k) * STEP;
      }
      p.push_back( f(s1) );
    }
    else if ( !(CGAL::is_negative(s0)) && !(CGAL::is_negative(s1)) )
    {
      FT tt = s0;
      int k = compute_k(tt, STEP);

      do
      {
        p.push_back( f(tt) );
        ++k;
        tt = FT(k * k) * STEP;
      }
      while ( CGAL::compare(tt, s0) == LARGER && CGAL::compare(tt, s1) == SMALLER );

      p.push_back( f(s1) );
    }
    else
    {
      FT tt = s1;
      int k = - compute_k(-tt, STEP);

      do
      {
        p.push_back( f(tt) );
        --k;
        tt = - FT(k * k) * STEP;
      }
      while ( CGAL::compare(tt, s0) == LARGER && CGAL::compare(tt, s1) == SMALLER );

      p.push_back( f(s0) );
    }
  }

  void generate_points(std::vector<Point_2>& p,
                       const FT STEP = FT(2)) const
  {
    return generate_points(p, STEP, t(p1), t(p2));
  }

  template< class Stream >
  void draw(Stream& W) const
  {
    std::vector< Point_2 > p;
    generate_points(p);

    for (unsigned int i = 0; i < p.size() - 1; i++) {
      W << Segment_2(p[i], p[i+1]);
    }
  }

  template< class K >
  void draw(CGAL::Qt::PainterOstream<K>& stream) const {
    stream.draw_parabola_segment(this->center(), this->line(), p1, p2);
  }
};



template< class Stream, class Gt >
inline
Stream& operator<<(Stream &s, const Parabola_segment_2<Gt> &P)
{
  P.draw(s);
  return s;
}

} //namespace CGAL

#endif // CGAL_PARABOLA_SEGMENT_2_H
