#ifndef CGAL_HYPERBOLA_2_H
#define CGAL_HYPERBOLA_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/enum.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/determinant.h>
#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

template < class Point, class Weight >
class Hyperbola_2
{
public:
  typedef CGAL::Weighted_point< Point, Weight >   Weighted_point;
  //  typedef typename R::RT          FT;
  typedef double                                  FT;
  typedef CGAL::Point_2< Cartesian<double>  >     Point_2;
  typedef CGAL::Segment_2< Cartesian< double > >  Segment_2;

protected:
  FT STEP;

  Point_2 f1, f2;
  FT r;
  Point_2 o;

  inline
  Point_2 lchain(const FT &t) const
  {
    std::vector< Point_2 > p = compute_points(t);
    if ( right(p[0]) )  return p[1];
    return p[0];
  }

  inline
  Point_2 rchain(const FT &t) const
  {
    std::vector< Point_2 > p = compute_points(t);
    if ( right(p[0]) )  return p[0];
    return p[1]; 
  }

  inline
  FT norm2(const Point_2& p) const
  {
    return (CGAL_NTS square(p.x()) + CGAL_NTS square(p.y()));
  }

  inline
  FT distance2(const Point_2& p1, const Point_2& p2) const
  {
    FT dx = p1.x()-p2.x();
    FT dy = p1.y()-p2.y();
    return (CGAL_NTS square(dx) + CGAL_NTS square(dy));
  }

  inline
  FT distance(const Point_2& p1, const Point_2& p2) const
  {
    return CGAL_NTS sqrt( distance2(p1, p2) );
  }

  void compute_origin()
  {
    FT dx = f2.x() - f1.x();
    FT dy = f2.y() - f1.y();
    FT a = CGAL_NTS sqrt(CGAL_NTS square(dx) + CGAL_NTS square(dy));
    FT t = (1 + r / a) / FT(2);

    o = Point_2(dx * t + f1.x(), dy * t + f1.y());
  }

  std::vector< Point_2 > compute_points(const FT &d) const {
      FT d1 = distance(o, f1) + d;
      FT d2 = distance(o, f2) + d;
      d1 *= d1;
      d2 *= d2;

      Point_2 df = Point_2(f2.x() - f1.x(), f2.y()-f1.y());

      std::vector< Point_2 > p;

      if ( CGAL_NTS is_negative(d) ) return p;

      if ( CGAL_NTS is_zero(df.x()) ) {
	FT y = (d1 - d2 + norm2(f2) - norm2(f1)) / (FT(2) * df.y());

	FT D = d1 - CGAL_NTS square(y - f1.y());

	D = CGAL_NTS abs(D);

	FT x1 = CGAL_NTS sqrt(D) + f1.x();
	FT x2 = -CGAL_NTS sqrt(D) + f1.x();

	p.push_back(Point_2(x1, y));
	p.push_back(Point_2(x2, y));

	return p;
      }

      FT gamma = (d1 - d2 + norm2(f2) - norm2(f1)) / (FT(2) * df.x());
      FT gamma1 = gamma - f1.x();
      FT beta = df.y() / df.x();

      FT a = FT(1) + CGAL_NTS square(beta);
      FT b = -FT(2) * (gamma1 * beta + f1.y());
      FT c = CGAL_NTS square(f1.y()) + CGAL_NTS square(gamma1) - d1;

      FT D = CGAL_NTS square(b) - FT(4) * a * c;

      D = CGAL_NTS abs(D);

      FT y1 = (-b + CGAL_NTS sqrt(D)) / (FT(2) * a);
      FT y2 = (-b - CGAL_NTS sqrt(D)) / (FT(2) * a);

      FT x1 = gamma - beta * y1;
      FT x2 = gamma - beta * y2;

      p.push_back(Point_2(x1, y1));
      p.push_back(Point_2(x2, y2));

      return p;
  }

  bool right(const Point_2& p) const
  {
    return CGAL_NTS
      is_negative( det3x3_by_formula< FT >(f1.x(), f1.y(), 1,
					   f2.x(), f2.y(), 1,
					   p.x(),  p.y(), 1) );
  }

  inline
  Point_2 midpoint(const Point_2& p1, const Point_2& p2) const
  {
    FT t1 = t(p1);
    FT t2 = t(p2);
    FT midt = (t1+t2)/2;
    return f(midt);
  }

  inline
  Point_2 f(FT t) const
  {
    if ( CGAL_NTS is_negative(t) )  return rchain(-t);
    return lchain(t);
  }

  inline
  FT t(const Point_2 &p) const
  {
    FT tt = distance(f1, p) - distance(f1, o);
    if ( right(p) )  return -tt;
    return tt;
  }

public:
  Hyperbola_2()  { STEP = FT(2); }

  Hyperbola_2(const Weighted_point &f1,	const Weighted_point &f2)
  {
    STEP = FT(2);
    this->r = CGAL_NTS to_double(f1.weight() - f2.weight());
    
    this->f1 = Point_2(CGAL_NTS to_double(f1.point().x()),
		       CGAL_NTS to_double(f1.point().y()));
    this->f2 = Point_2(CGAL_NTS to_double(f2.point().x()),
		       CGAL_NTS to_double(f2.point().y()));

    compute_origin();
  }

  Oriented_side
  side_of_hyperbola(const Point &p) const
  {
    Point_2 q(CGAL_NTS to_double(p.x()),
	      CGAL_NTS to_double(p.y()));
    double dist = distance(q, f1) - distance(q, f2) - r;
    if ( dist < 0 )  return ON_NEGATIVE_SIDE;
    if ( dist > 0 )  return ON_POSITIVE_SIDE;
    return ON_ORIENTED_BOUNDARY;
  }


#if defined CGAL_QT_WIDGET_H
  void draw_qt(Qt_widget& W) const
    {
      std::vector< Point_2 > p;
      std::vector< Point_2 > pleft, pright;

      pleft.push_back(o);
      pright.push_back(o);

      double width = W.x_max() - W.x_min();
      double height = W.y_max() - W.y_min();

      FT STEP;
      if ( width < height ) {
	STEP = width / 100.0;
      } else {
	STEP = height / 100.0;
      }
      //    double mind = distance(o, f1) - r1;
      for (int i = 1; i <= 100; i++) {
	p = compute_points(i * i * STEP);
	
	if ( p.size() > 0 ) {
	  if ( right(p[0]) ) {
	    pright.push_back(p[0]);
	    pleft.push_back(p[1]);
	  } else {
	    pright.push_back(p[1]);
	    pleft.push_back(p[0]);
	  }
	}
      }

      for (unsigned int i = 0; i < pleft.size() - 1; i++) {
	W << Segment_2(pleft[i], pleft[i+1]);
      }
      
      for (unsigned int i = 0; i < pright.size() - 1; i++) {
	W << Segment_2(pright[i], pright[i+1]);
      }
    }
#endif

  template< class Stream >
  void draw(Stream &W) const
  {
    std::vector< Point_2 > p;
    std::vector< Point_2 > pleft, pright;

    pleft.push_back(o);
    pright.push_back(o);

    //    double mind = distance(o, f1) - r1;
    for (int i = 1; i <= 100; i++) {
      p = compute_points(i * i * STEP);

      if ( p.size() > 0 ) {
	if ( right(p[0]) ) {
	  pright.push_back(p[0]);
	  pleft.push_back(p[1]);
	} else {
	  pright.push_back(p[1]);
	  pleft.push_back(p[0]);
	}
      }
    }

    for (unsigned int i = 0; i < pleft.size() - 1; i++) {
      W << Segment_2(pleft[i], pleft[i+1]);
    }

    for (unsigned int i = 0; i < pright.size() - 1; i++) {
      W << Segment_2(pright[i], pright[i+1]);
    }
  }

};

template< class Stream, class Point, class Weight >
inline
Stream& operator<<(Stream& s,
		   const Hyperbola_2< Point, Weight > &H)
{
  H.draw(s);
  return s;
}

#if defined CGAL_QT_WIDGET_H
template< class Point, class Weight >
inline
Qt_widget& operator<<(Qt_widget& s,
		      const Hyperbola_2< Point, Weight > &H)
{
  H.draw_qt(s);
  return s;
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_HYPERBOLA_2_H
