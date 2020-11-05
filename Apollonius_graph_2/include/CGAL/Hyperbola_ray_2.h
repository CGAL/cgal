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



#ifndef CGAL_HYPERBOLA_RAY_2_H
#define CGAL_HYPERBOLA_RAY_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Hyperbola_segment_2.h>

namespace CGAL {

template < class Gt >
class Hyperbola_ray_2 : public Hyperbola_segment_2< Gt >
{
public:
  typedef Sign                             Hyperbola_direction;
  typedef CGAL::Hyperbola_segment_2<Gt>    Base;
  typedef typename Base::Site_2            Site_2;
  typedef typename Base::Point_2           Point_2;
  typedef typename Base::Segment_2         Segment_2;
  typedef typename Gt::Ray_2               Ray_2;
  typedef typename Base::FT                FT;
  //  typedef typename R::RT         FT;
  //  typedef double                                 FT;
  //  typedef CGAL::Point_2< Cartesian<double> >     Point_2;
  //  typedef CGAL::Segment_2< Cartesian<double> >   Segment_2;
  //  typedef CGAL::Ray_2< Cartesian<double> >       Ray_2;


  using Base::t;
  using Base::f;

protected:
#if defined(__POWERPC__) && \
  defined(__GNUC__) && (__GNUC__ == 3 ) && (__GNUC_MINOR__ == 4)
  // hack to avoid nasty warning for G++ 3.4 on Darwin
  static FT OFFSET()
  {
    return FT(10000);
  }
#else
  static const FT& OFFSET()
  {
    static const FT offset_(10000);
    return offset_;
  }
#endif

  template< class Stream >
  inline
  void draw_ray(Stream &W) const
  {
    W << Ray_2(this->p1, this->p2);
  }

  Site_2 _f1, _f2;
  Point_2 _p;
  Hyperbola_direction _dir;

public:
  Hyperbola_ray_2() : Hyperbola_segment_2< Gt >() {}


  Hyperbola_ray_2(const Site_2 &f1, const Site_2 &f2,
                  const Point_2 &p,
                  const Hyperbola_direction& direction) :
    Hyperbola_segment_2< Gt >(f1, f2, p, p),
    _f1(f1), _f2(f2), _p(p), _dir(direction)
  {
    FT t1 = t(this->p1);
    if ( direction == POSITIVE ) {
      this->p2 = f(t1 + this->STEP * OFFSET());
    } else {
      this->p2 = f(t1 - this->STEP * OFFSET());
    }
  }


  template<class QTWIDGET>
  void draw_qt(QTWIDGET& s)
  {
    if ( CGAL::is_zero(this->r) ) {
      draw_ray(s);
      return;
    }

    double width = s.x_max() - s.x_min();
    double height = s.y_max() - s.y_min();

    if ( width > height ) {
      this->STEP = height / 100.0;
    } else {
      this->STEP = width / 100.0;
    }

    FT t1 = t(this->p1);
    if ( _dir == POSITIVE ) {
      this->p2 = f(t1 + this->STEP * OFFSET());
    } else {
      this->p2 = f(t1 - this->STEP * OFFSET());
    }

    Hyperbola_segment_2< Gt >::draw(s);
  }

  template< class Stream >
  inline
  void draw(Stream& s) const
  {
    if ( CGAL::is_zero(this->r) ) {
      draw_ray(s);
      return;
    }

    Hyperbola_segment_2< Gt >::draw(s);
  }

};



template< class Stream, class Gt >
inline
Stream& operator<<(Stream &s, const Hyperbola_ray_2<Gt> &H)
{
  H.draw(s);
  return s;
}

} //namespace CGAL

#endif // CGAL_HYPERBOLA_RAY_2_H
