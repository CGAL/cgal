// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_PARABOLA_SEGMENT_2_H
#define CGAL_PARABOLA_SEGMENT_2_H

#include <CGAL/Parabola_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt >
class Parabola_segment_2 : public Parabola_2< Gt >
{
  typedef CGAL::Parabola_2<Gt>            Base;
  typedef typename Base::Site_2           Site_2;
  typedef typename Base::FT               FT;
  typedef typename Base::Point_2          Point_2;
  typedef typename Base::Segment_2        Segment_2;
  typedef typename Base::Line_2           Line_2;

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

  void generate_points(std::vector<Point_2>& p) const
  {
    FT s[2];

    s[0] = t(p1);
    s[1] = t(p2);

    if (CGAL::compare(s[0], s[1]) == LARGER) {
#if defined(__GNUC__) && (__GNUC__ < 3)
      FT tmp = s[0];
      s[0] = s[1];
      s[1] = tmp;
#else
      std::swap< FT >(s[0], s[1]);
#endif
    }

    p.clear();

    if ( !(CGAL::is_positive(s[0])) &&
	 !(CGAL::is_negative(s[1])) ) {
      FT tt;
      int k;

      p.push_back( this->o );
      k = 1;
      tt = -this->STEP();
      while ( CGAL::compare(tt, s[0]) == LARGER ) {
	p.insert( p.begin(), f(tt) );
	k--;
	tt = -FT(k * k) * this->STEP();
      }
      p.insert( p.begin(), f(s[0]) );

      k = 1;
      tt = this->STEP();
      while ( CGAL::compare(tt, s[1]) == SMALLER ) {
	p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * this->STEP();
      }
      p.push_back( f(s[1]) );
    } else if ( !(CGAL::is_negative(s[0])) &&
		!(CGAL::is_negative(s[1])) ) {
      FT tt;
      int k;


      p.push_back( f(s[0]) );

      tt = s[0];
      k = int(CGAL::to_double(CGAL::sqrt(tt / this->STEP())));

      while ( CGAL::compare(tt, s[1]) == SMALLER ) {
	if ( CGAL::compare(tt, s[0]) != SMALLER )
	  p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * this->STEP();
      }
      p.push_back( f(s[1]) );
    } else {
      FT tt;
      int k;

      p.push_back( f(s[1]) );

      tt = s[1];
      k = int(CGAL::to_double(-CGAL::sqrt(-tt / this->STEP())));

      while ( CGAL::compare(tt, s[0]) == LARGER ) {
	if ( CGAL::compare(tt, s[1]) != LARGER )
	  p.push_back( f(tt) );
	k--;
	tt = -FT(k * k) * this->STEP();
      }
      p.push_back( f(s[0]) );
    }
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
};



template< class Stream, class Gt >
inline
Stream& operator<<(Stream &s, const Parabola_segment_2<Gt> &P)
{
  P.draw(s);
  return s;
}

CGAL_END_NAMESPACE

#endif // CGAL_PARABOLA_SEGMENT_2_H
