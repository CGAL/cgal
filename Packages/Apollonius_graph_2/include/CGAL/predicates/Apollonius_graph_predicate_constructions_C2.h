#ifndef CGAL_APOLLONIUS_GRAPH_PREDICATE_CONSTRUCTIONS_C2_H
#define CGAL_APOLLONIUS_GRAPH_PREDICATE_CONSTRUCTIONS_C2_H

#include <CGAL/Weighted_point.h>

CGAL_BEGIN_NAMESPACE

template< class Point, class Weight >
class Inverted_weighted_point
  : public Weighted_point<Point, Weight>
{
public:
  typedef Weighted_point<Point, Weight>    Weighted_point;
  typedef typename Point::R::FT            FT;
private:
  FT   _p;
public:
  Inverted_weighted_point(const Weighted_point& wp, const FT& p)
    : Weighted_point(wp), _p(p) {}

  inline FT p() const { return _p; }
};


template< class Point, class Weight >
class Weighted_point_inverter
{
public:
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef typename Point::R::FT                   FT;
private:
  Weighted_point _pole;
public:
  Weighted_point_inverter(const Weighted_point& pole)
    : _pole(pole) {}

  Inverted_weighted_point operator()(const Weighted_point& wp)
    {
      FT xs = wp.x() - _pole.x();
      FT ys = wp.y() - _pole.y();
      FT ws = wp.weight() - _pole.weight();
      FT ps = CGAL_NTS square(xs) + CGAL_NTS square(ys)
	- CGAL_NTS square(ws);

      return
	Inverted_weighted_point(Weighted_point(Point(xs, ys), ws), ps);
    }

  Weighted_point pole() const { return _pole; }
};


template< class Point, class Weight >
class Voronoi_radius
{
  // this class stores the coefficients for the tritangent circle
  // radius equation. In particular we have:
  //             a x^2 - 2 b x + c = 0;
  // x here represents the inverse of the radius
public:
  typedef Weighted_point<Point, Weight>            Weighted_point;
  typedef typename Point::R::FT                    FT;
  typedef Inverted_weighted_point<Point, Weight>   Inverted_weighted_point;
private:
  FT _a, _b, _c;
  FT _c2, _delta;
  FT _dxp, _dyp, _dwp;

  Voronoi_radius(FT a, FT b, FT c, FT c2, FT delta,
		 FT dxp, FT dyp, FT dwp)
    : _a(a), _b(b), _c(c), _c2(c2), _delta(delta), _dxp(dxp),
      _dyp(dyp), _dwp(dwp) {}

public:
  Voronoi_radius(const Inverted_weighted_point& u1,
		 const Inverted_weighted_point& u2)
    {
      FT dxp = det2x2_by_formula(u1.x(), u1.p(), u2.x(), u2.p());
      FT dyp = det2x2_by_formula(u1.y(), u1.p(), u2.y(), u2.p());
      FT dwp = det2x2_by_formula(u1.weight(), u1.p(), u2.weight(), u2.p());
      FT dxy = det2x2_by_formula(u1.x(), u1.y(), u2.x(), u2.y());
      FT dxw = det2x2_by_formula(u1.x(), u1.weight(), u2.x(), u2.weight());
      FT dyw = det2x2_by_formula(u1.y(), u1.weight(), u2.y(), u2.weight());

      _a = CGAL_NTS square(dxp) + CGAL_NTS square(dyp);
      _b = dxp * dxw + dyp * dyw;
      _c = CGAL_NTS square(dxw) + CGAL_NTS square(dyw) - CGAL_NTS square(dxy);
      _c2 = dxy;
      _delta = _a - CGAL_NTS square(dwp);
      _dxp = dxp;
      _dyp = dyp;
      _dwp = dwp;
    }

  inline FT a() const { return _a; }
  inline FT b() const { return _b; }
  inline FT c() const { return _c; }

  inline FT c1() const { return _b; }
  inline FT c2() const { return _c2; }
  inline FT delta() const { return _delta; }
  inline FT d() const { return _a; }

  inline FT dxp() const { return _dxp; }
  inline FT dyp() const { return _dyp; }
  inline FT dwp() const { return _dwp; }

  inline bool is_first_root() const { return CGAL_NTS is_negative(_c2); }

  Voronoi_radius get_symmetric()
    {
      return Voronoi_radius(_a, _b, _c, -_c2, _delta, -_dxp, -_dyp, -_dwp);
    }
};


template< class Point, class Weight >
class Bitangent_line
{
  // this class computes and stores the data for the left bitangent
  // line of the weighted points p1, p2 oriented from p1 to p2
  // or the left bitangent line of the inverted weighted point u1 and
  // u2, oriented from u1 to u2
public:
  typedef Weighted_point<Point, Weight>            Weighted_point;
  typedef Inverted_weighted_point<Point, Weight>   Inverted_weighted_point;
  typedef typename Point::R::FT                    FT;
protected:
  FT _a1, _a2;
  FT _b1, _b2;
  FT _c1, _c2;
  FT _delta;
  FT _d;
  FT _dw;
  FT _dxw, _dyw;

  Bitangent_line(FT a1, FT a2, FT b1, FT b2, FT c1, FT c2,
		 FT delta, FT d, FT dw, FT dxw, FT dyw)
    : _a1(a1), _a2(a2), _b1(b1), _b2(b2), _c1(c1), _c2(c2),
      _delta(delta), _d(d), _dw(dw),_dxw(dxw), _dyw(dyw) {}

  inline void
  store(FT dx, FT dy, FT dw)
    {
      _dw = dw; 
      _a1 = dx * dw;
      _a2 = dy;
      _b1 = dy * dw;
      _b2 = -dx;
    }

  inline void
  store(FT dx, FT dy, FT dw, FT dxy, FT dxw, FT dyw)
    {
      store(dx, dy, dw);
      _c1 = dx * dxw + dy * dyw;
      _c2 = dxy;
      _d = CGAL_NTS square(dx) + CGAL_NTS square(dy);
      _delta = _d - CGAL_NTS square(dw);
      _dxw = dxw;
      _dyw = dyw;
    }

public:
  Bitangent_line(const Weighted_point& p1, const Weighted_point& p2)
    {
      FT dx = p1.x() - p2.x();
      FT dy = p1.y() - p2.y();
      FT dw = p1.weight() - p2.weight();
      FT dxy = det2x2_by_formula(p1.x(), p1.y(), p2.x(), p2.y());
      FT dxw = det2x2_by_formula(p1.x(), p1.weight(), p2.x(), p2.weight());
      FT dyw = det2x2_by_formula(p1.y(), p1.weight(), p2.y(), p2.weight());

      store(dx, dy, dw, dxy, dxw, dyw);
    }


  Bitangent_line(const Inverted_weighted_point& u1,
		 const Inverted_weighted_point& u2)
    {
      FT dxp = det2x2_by_formula(u1.x(), u1.p(), u2.x(), u2.p());
      FT dyp = det2x2_by_formula(u1.y(), u1.p(), u2.y(), u2.p());
      FT dwp = det2x2_by_formula(u1.weight(), u1.p(), u2.weight(), u2.p());
      FT dxy = det2x2_by_formula(u1.x(), u1.y(), u2.x(), u2.y());
      FT dxw = det2x2_by_formula(u1.x(), u1.weight(), u2.x(), u2.weight());
      FT dyw = det2x2_by_formula(u1.y(), u1.weight(), u2.y(), u2.weight());

      store(dxp, dyp, dwp, dxy, dxw, dyw);
    }

  Bitangent_line get_symmetric() const
    {
      return
	Bitangent_line(_a1, -_a2, _b1, -_b2, _c1, -_c2, _delta, _d,
		       -_dw, -_dxw, -_dyw);
    }

  Bitangent_line get_rot90() const
    {
      return
	Bitangent_line(-_b1, -_b2, _a1, _a2, _c1, _c2, _delta, _d,
		       _dw, -_dyw, _dxw);
    }

  Bitangent_line perpendicular(const Point& p) const
    {
      // THIS DOES NOT KEEP TRACK OF THE ADDITIONALLY STORED
      // QUANTITIES; THIS IS INEVITABLE IN ANY CASE SINCE GIVEN p WE
      // CANNOT ANY LONGER HOPE TO KEEP TRACK OF THOSE
      Bitangent_line rotated = get_rot90();
      rotated._c1 = _b1 * p.x() - _a1 * p.y();
      rotated._c2 = _b2 * p.x() - _a2 * p.y();

      return rotated;
    }

  Bitangent_line perpendicular(const Inverted_weighted_point& u) const
    {
      // THIS DOES NOT KEEP TRACK OF THE ADDITIONALLY STORED
      // QUANTITIES; THIS IS INEVITABLE IN ANY CASE SINCE GIVEN p WE
      // CANNOT ANY LONGER HOPE TO KEEP TRACK OF THOSE
      Bitangent_line rotated = get_rot90();
      rotated._c1 = (_b1 * u.x() - _a1 * u.y()) * u.p();
      rotated._c2 = (_b2 * u.x() - _a2 * u.y()) * u.p();

      return rotated;
    }

  inline FT a1() const { return _a1; }
  inline FT a2() const { return _a2; }
  inline FT b1() const { return _b1; }
  inline FT b2() const { return _b2; }
  inline FT c1() const { return _c1; }
  inline FT c2() const { return _c2; }
  inline FT delta() const { return _delta; }
  inline FT d() const { return _d; }

  inline FT dx() const { return -_b2; }
  inline FT dy() const { return _a2; }
  inline FT dw() const { return _dw; }
  inline FT dxw() const { return _dxw; }
  inline FT dyw() const { return _dyw; }
};


template<class Point, class Weight>
class Voronoi_circle : public Bitangent_line<Point, Weight>
{
public:
  //  typedef Weighted_point<Point, Weight>            Weighted_point;
  typedef CGAL::Inverted_weighted_point<Point,Weight>  Inverted_weighted_point;
  typedef Bitangent_line<Point, Weight>            Bitangent_line;
  //  typedef Bitangent_line::Inverted_weighted_point  Inverted_weighted_point;
  typedef Voronoi_radius<Point, Weight>            Voronoi_radius;
  //  typedef typename Point::R::FT                    FT;
protected:
  FT _gamma;

  inline
  void compute_gamma()
    {
      _gamma = CGAL_NTS square(_dxw) + CGAL_NTS square(_dyw)
	- CGAL_NTS square(_c2);
    }

public:
  Voronoi_circle(const Voronoi_radius& vr)
    : Bitangent_line(FT(0), FT(0), FT(0), FT(0), vr.b(), vr.c2(),
		     vr.delta(), vr.d(), FT(0), FT(0), FT(0)), _gamma(vr.c())
    {
      store(vr.dxp(), vr.dyp(), vr.dwp());
    }

  Voronoi_circle(const Bitangent_line& bl)
    : Bitangent_line(bl.a1(), bl.a2(), bl.b1(), bl.b2(), bl.c1(), bl.c2(),
		     bl.delta(), bl.d(), bl.dw(), bl.dxw(), bl.dyw())
    {
      compute_gamma();
    }
		
  inline FT alpha() const { return _d; }
  inline FT beta() const { return  _c1; }
  inline FT gamma() const { return _gamma; }

  inline bool is_first_root() const { return CGAL_NTS is_negative(_c2); }

  FT compute_P4(const Inverted_weighted_point& u1,
		const Inverted_weighted_point& u2,
		const Inverted_weighted_point& u3) const
    {
      FT dx1 = det2x2_by_formula(u2.x(), u2.p(), u1.x(), u1.p());
      FT dy1 = det2x2_by_formula(u2.y(), u2.p(), u1.y(), u1.p());
      FT dw1 = det2x2_by_formula(u2.weight(), u2.p(), u1.weight(), u1.p());

      FT dx3 = det2x2_by_formula(u3.x(), u3.p(), u2.x(), u2.p());
      FT dy3 = det2x2_by_formula(u3.y(), u3.p(), u2.y(), u2.p());
      FT dw3 = det2x2_by_formula(u3.weight(), u3.p(), u2.weight(), u2.p());

      FT u2Pv2 = CGAL_NTS square(u2.x()) + CGAL_NTS square(u2.y());
      FT u2Mv2 = CGAL_NTS square(u2.x()) - CGAL_NTS square(u2.y());
      FT u2v2 = FT(2) * u2.x() * u2.y();


      FT vvMuu = dy1 * dy3 - dx1 * dx3;
      FT vuPuv = dy1 * dx3 + dx1 * dy3;

      FT dx2Pdy2_1 = CGAL_NTS square(dx1) + CGAL_NTS square(dy1);
      FT dx2Pdy2_3 = CGAL_NTS square(dx3) + CGAL_NTS square(dy3);

      FT fr1_sq = CGAL_NTS square(dw1) * dx2Pdy2_3;
      FT fr3_sq = CGAL_NTS square(dw3) * dx2Pdy2_1;

      FT f1 = (fr1_sq + fr3_sq) * CGAL_NTS square(u2Pv2);

      FT f2 = FT(2) * dw1 * dw3 * u2Pv2 * (u2Mv2 * vvMuu - u2v2 * vuPuv );
      FT f3 = CGAL_NTS square(u2Mv2 * vuPuv + u2v2 * vvMuu);

      FT F = f1 + f2 - f3;


      FT uuPvv = dy1 * dy3 + dx1 * dx3;
      FT vuMuv = dy1 * dx3 - dx1 * dy3;

      FT G = fr1_sq + fr3_sq - FT(2) * dw1 * dw3 * uuPvv
	- CGAL_NTS square(vuMuv);

      return (F * G);
    }
};



CGAL_END_NAMESPACE

#endif  // CGAL_APOLLONIUS_GRAPH_PREDICATE_CONSTRUCTIONS_C2_H
