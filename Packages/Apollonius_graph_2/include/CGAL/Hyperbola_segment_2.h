#ifndef CGAL_HYPERBOLA_SEGMENT_2_H
#define CGAL_HYPERBOLA_SEGMENT_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Hyperbola_2.h>

CGAL_BEGIN_NAMESPACE

template < class Point, class Weight >
class Hyperbola_segment_2 : public Hyperbola_2< Point, Weight >
{
public:
  typedef CGAL::Weighted_point< Point, Weight >   Weighted_point;
  //  typedef typename R::RT          FT;
  typedef double                                  FT;
  typedef CGAL::Point_2< Cartesian<double> >      Point_2;
  typedef CGAL::Segment_2< Cartesian<double> >    Segment_2;

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
    return Hyperbola_2< Point, Weight >::midpoint(p1, p2);
  }

public:
  Hyperbola_segment_2() : Hyperbola_2< Point, Weight >() {}

  Hyperbola_segment_2(const Weighted_point &f1,
		      const Weighted_point &f2,
		      const Point &p1,
		      const Point &p2) :
    Hyperbola_2< Point, Weight >(f1, f2)
  {
    this->p1 = Point_2(CGAL_NTS to_double(p1.x()),
		       CGAL_NTS to_double(p1.y()));
    this->p2 = Point_2(CGAL_NTS to_double(p2.x()),
		       CGAL_NTS to_double(p2.y()));
  }

  template< class Stream >
  void draw(Stream &W) const
  {
    if ( CGAL_NTS is_zero(r) ) {
      draw_line(W);
      return;
    }

    vector< Point_2 > p;

    //    FT STEP = W.width() / 100.0;

    FT s[2];

    s[0] = t(p1);
    s[1] = t(p2);

    if (CGAL_NTS compare(s[0], s[1]) == LARGER)
      swap< FT >(s[0], s[1]);

    p.clear();

    if ( !(CGAL_NTS is_positive(s[0])) &&
	 !(CGAL_NTS is_negative(s[1])) ) {
      FT tt;
      int k;

      p.push_back( o );
      k = 1;
      tt = FT(-STEP);
      while ( CGAL_NTS compare(tt, s[0]) == LARGER ) {
	p.insert( p.begin(), f(tt) );
	k--;
	tt = -FT(k * k) * STEP;
      }
      p.insert( p.begin(), f(s[0]) );

      k = 1;
      tt = FT(STEP);
      while ( CGAL_NTS compare(tt, s[1]) == SMALLER ) {
	p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * STEP;
      }
      p.push_back( f(s[1]) );
    } else if ( !(CGAL_NTS is_negative(s[0])) &&
		!(CGAL_NTS is_negative(s[1])) ) {
      FT tt;
      int k;


      p.push_back( f(s[0]) );

      tt = s[0];
      k = int(CGAL_NTS to_double(CGAL_NTS sqrt(tt / STEP)));

      while ( CGAL_NTS compare(tt, s[1]) == SMALLER ) {
	if ( CGAL_NTS compare(tt, s[0]) != SMALLER )
	  p.push_back( f(tt) );
	k++;
	tt = FT(k * k) * STEP;
      }
      p.push_back( f(s[1]) );
    } else {
      FT tt;
      int k;

      p.push_back( f(s[1]) );

      tt = s[1];
      k = int(CGAL_NTS to_double(-CGAL_NTS sqrt(-tt / STEP)));

      while ( CGAL_NTS compare(tt, s[0]) == LARGER ) {
	if ( CGAL_NTS compare(tt, s[1]) != LARGER )
	  p.push_back( f(tt) );
	k--;
	tt = -FT(k * k) * STEP;
      }
      p.push_back( f(s[0]) );
    }

    for (unsigned int i = 0; i < p.size() - 1; i++) {
      W << Segment_2(p[i], p[i+1]);
    }
  }
};


template< class Stream, class Point, class Weight >
inline
Stream& operator<<(Stream &s,
		   const Hyperbola_segment_2< Point, Weight > &H)
{
  H.draw(s);
  return s;
}

CGAL_END_NAMESPACE

#endif // CGAL_HYPERBOLA_SEGMENT_2_H
