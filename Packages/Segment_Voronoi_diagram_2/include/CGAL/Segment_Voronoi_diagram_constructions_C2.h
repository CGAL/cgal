#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_CONSTRUCTIONS_C2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_CONSTRUCTIONS_C2_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>

#include <CGAL/predicates/Square_root_1.h>
#include <CGAL/predicates/Square_root_2.h>


#include <CGAL/Parabola_2.h>
#include <CGAL/Parabola_segment_2.h>

CGAL_BEGIN_NAMESPACE


template< class K >
class Svd_voronoi_vertex_2
{
public:
  typedef enum {PPP = 0, PPS, PSS, SSS} vertex_t;
  struct PPP_Type {};
  struct PPS_Type {};
  struct PSS_Type {};
  struct SSS_Type {};

  typedef typename K::RT         RT;
  typedef typename K::FT         FT;
  typedef typename K::Point_2    Point;
  typedef typename K::Circle_2   Circle;
  typedef typename K::Segment_2  Segment;
  typedef typename K::Line_2     Line;
  typedef typename K::Site_2     Site;

  typedef CGAL::Square_root_1<RT>       Sqrt_1;
  typedef CGAL::Square_root_2<RT>       Sqrt_2;
  typedef CGAL::Square_root_2<Sqrt_1>   Sqrt_3;


  struct Homogeneous_point_2
  {
    RT hx, hy, hw;

    Homogeneous_point_2() : hx(0), hy(0), hw(0) {}

    FT x() const { return hx / hw; }
    FT y() const { return hy / hw; }
  };

public:
  static
  void compute_supporting_line(const Segment& s, RT& a, RT& b, RT& c)
  {
    a = s.source().y() - s.target().y();
    b = s.target().x() - s.source().x();
    c = s.source().x() * s.target().y() - s.target().x() * s.source().y();
  }

  static
  Homogeneous_point_2
  compute_projection(const Line& l, const Point& p)
  {
    Homogeneous_point_2 hp;

    RT ab = l.a() * l.b();

    hp.hx = CGAL_NTS square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    hp.hy = CGAL_NTS square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    hp.hw = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

    return hp;
  }


private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Point& p, const Point& q, const Point& r)
  {
    v_type = PPP;

    RT np = CGAL_NTS square(p.x()) + CGAL_NTS square(p.y());
    RT nq = CGAL_NTS square(q.x()) + CGAL_NTS square(q.y());
    RT nr = CGAL_NTS square(r.x()) + CGAL_NTS square(r.y());

    ux = np * (q.y() - r.y()) + nq * (r.y() - p.y()) + nr * (p.y() - q.y());
    uy = -(np * (q.x() - r.x()) + nq * (r.x() - p.x()) + nr * (p.x() - q.x()));
    uz = RT(2) * ( (q.x() * r.y() - r.x() * q.y()) +
		   (r.x() * p.y() - p.x() * r.y()) +
		   (p.x() * q.y() - q.x() * p.y()) );
  }

  //--------------------------------------------------------------------------

  void
  compute_pss(const Point& p, const Segment& q, const Segment& r)
  {
    v_type = PSS;

    _s = ZERO;
    _r = ZERO;

    bool pq = (p == q.source() || p == q.target());
    bool pr = (p == r.source() || p == r.target());

    if ( pq && pr ) {
      J = RT(0);
      I = RT(0);
      c1c2 = RT(0);
      D1D2 = RT(0);
      a1a2 = RT(-1);
      b1b2 = RT(0);
      return;
    }


    RT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(q, a1, b1, c1);
    compute_supporting_line(r, a2, b2, c2);

    RT c1_ = a1 * p.x() + b1 * p.y() + c1;
    RT c2_ = a2 * p.x() + b2 * p.y() + c2;

    if ( p == q.source() || p == q.target() ) {
      c1_ = RT(0);
    }
    if ( p == r.source() || p == r.target() ) {
      c2_ = RT(0);
    }

    Sign sgn_c1_ = CGAL_NTS sign(c1_);
    Sign sgn_c2_ = CGAL_NTS sign(c2_);

    if ( sgn_c1_ == NEGATIVE ) {
      a1 = -a1;  b1 = -b1;  c1_ = -c1_;
    } else if ( sgn_c1_ == ZERO ) {
      CGAL_assertion( p == q.source() || p == q.target() );
      if ( p == q.target() ) {
	a1 = -a1;  b1 = -b1;  c1_ = -c1_;
      }
    }

    if ( sgn_c2_ == NEGATIVE ) {
      a2 = -a2;  b2 = -b2;  c2_ = -c2_;
    } else if ( sgn_c2_ == ZERO ) {
      CGAL_assertion( p == r.source() || p == r.target() );
      if ( p == r.source() ) {
	a2 = -a2;  b2 = -b2;  c2_ = -c2_;
      }
    }

    if ( pq ) {
      J = a1 * c2_;
      I = b1 * c2_;
      a1a2 = a1 * a2;
      b1b2 = b1 * b2;
      c1c2 = RT(0);

      RT n1 = CGAL_NTS square(a1) + CGAL_NTS square(b1);
      RT n2 = CGAL_NTS square(a2) + CGAL_NTS square(b2);

      D1D2 = n1 * n2;
    } else if ( pr ) {
      J = a2 * c1_;
      I = b2 * c1_;
      a1a2 = a1 * a2;
      b1b2 = b1 * b2;
      c1c2 = RT(0);

      RT n1 = CGAL_NTS square(a1) + CGAL_NTS square(b1);
      RT n2 = CGAL_NTS square(a2) + CGAL_NTS square(b2);

      D1D2 = n1 * n2;
    } else {
      Line lq(a1, b1, c1_);
      Line lr(a2, b2, c2_);
      compute_pll(p, lq, lr);
    }
  }


  void
  compute_pll(const Point& p, const Line& lq, const Line& lr)
  {
    RT a1 = lq.a(), b1 = lq.b(), c1_ = lq.c();
    RT a2 = lr.a(), b2 = lr.b(), c2_ = lr.c();

    CGAL_precondition( c1_ >= 0 );
    CGAL_precondition( c2_ >= 0 );

    RT n1 = CGAL_NTS square(a1) + CGAL_NTS square(b1);
    RT n2 = CGAL_NTS square(a2) + CGAL_NTS square(b2);

    I = b1 * c2_ + b2 * c1_;
    J = a1 * c2_ + a2 * c1_;

    c1c2 = RT(2) * c1_ * c2_;
    a1a2 = a1 * a2;
    b1b2 = b1 * b2;

    D1D2 = n1 * n2;

    Sign s1, s2;

    s1 = CGAL_NTS sign(b1);
    s2 = CGAL_NTS sign(-b2);
    if ( s1 == ZERO ) {
      _s = s2;
    } else if ( s2 == ZERO ) {
      _s = s1;
    } else if ( s1 == s2 ) {
      _s = s1;
    } else {
      _s = CGAL_NTS sign(b1);

      RT e = b1 * b1 * n2 - b2 * b2 * n1;
      _s = Sign(s1 * CGAL_NTS sign(e));
#if 0
      RT e1 = a2 * b1 - a1 * b2;
      RT e2 = a1 * b2 + a2 * b1;

      _s = Sign(_s * CGAL_NTS sign(e1) * CGAL_NTS sign(e2));
#endif
    }

    s1 = CGAL_NTS sign(a1);
    s2 = CGAL_NTS sign(-a2);
    if ( s1 == ZERO ) {
      _r = s2;
    } else if ( s2 == ZERO ) {
      _r = s1;
    } else if ( s1 == s2 ) {
      _r = s1;
    } else {
      _r = s1;

      RT e = a1 * a1 * n2 - a2 * a2 * n1;
      _r = Sign(s1 * CGAL_NTS sign(e));
#if 0
      RT e1 = a1 * b2 - a2 * b1;
      RT e2 = a1 * b2 + a2 * b1;

      _r = Sign(_r * CGAL_NTS sign(e1) * CGAL_NTS sign(e2));
#endif
    }
  }


  //--------------------------------------------------------------------------


  void
  compute_pps(const Point& p, const Point& q, const Segment& r)
  {
    v_type = PPS;

    RT a, b, c;
    compute_supporting_line(r, a, b, c);

    RT c_ = a * p.x() + b * p.y() + c;
    RT cq_ = a * q.x() + b * q.y() + c;


    if ( p == r.source() || p == r.target() ) {
      c_ = RT(0);
    }
    if ( q == r.source() || q == r.target() ) {
      cq_ = RT(0);
    }

    Sign s = CGAL_NTS sign(c_);

    if ( s == NEGATIVE ) {
      a = -a;  b = -b;  c = -c;  c_ = -c_;  cq_ = -cq_;
    } else if ( s == ZERO ) {
      Sign s1 = CGAL_NTS sign(cq_);

      CGAL_assertion( s1 != ZERO );
      if ( s1 == NEGATIVE ) {
	a = -a;  b = -b;  c = -c;  c_ = -c_;  cq_ = -cq_;
      }
    }

    RT nl = CGAL_NTS square(a) + CGAL_NTS square(b);

    x_ = q.x() - p.x();
    y_ = q.y() - p.y();
    RT n_ = CGAL_NTS square(x_) + CGAL_NTS square(y_);


    Comparison_result res = CGAL_NTS compare( c_, cq_ );

    if ( res == EQUAL ) {
      RT e1 = CGAL_NTS square(c_);
      J = nl * (a * n_ + RT(4) * c_ * x_) - RT(4) * a * e1;
      I = nl * (b * n_ + RT(4) * c_ * y_) - RT(4) * b * e1;
      X = RT(8) * nl * c_;

      // have to reset x_ and y_ so that the computation of the
      // Voronoi vertex can be done in a unified way for the two cases.
      x_ = y_ = RT(0);
      S = RT(0);
      return;
    }


    RT e1 = a * x_ + b * y_;
    RT e2 = b * x_ - a * y_;
    RT e3 = n_ * e1;
    RT e4 = RT(2) * c_ * e2;

    X = RT(2) * CGAL_NTS square(e1);
    I = b * e3 + x_ * e4;
    J = a * e3 - y_ * e4;
    S = n_ * nl * c_ * cq_;
  }


  //--------------------------------------------------------------------------

  bool is_on_positive_halfspace(const Line& l, const Segment& s) const
  {
    Oriented_side os1, os2;

    os1 = l.oriented_side(s.source());
    os2 = l.oriented_side(s.target());

    return ( (os1 == ON_POSITIVE_SIDE && os2 != ON_NEGATIVE_SIDE) ||
	     (os1 != ON_NEGATIVE_SIDE && os2 == ON_POSITIVE_SIDE) );
  }

  bool
  is_consistent(const Segment& p, const Segment& q, const Segment& r,
		RT a[], RT b[], RT c[]) const
  {
    Line l[3] = {Line(a[0], b[0], c[0]), Line(a[1], b[1], c[1]), 
		 Line(a[2], b[2], c[2])};

    int num_oriented(0);

    if ( is_on_positive_halfspace(l[0], q) ||
	 is_on_positive_halfspace(l[0], r) ) {
      num_oriented++;
    }

    if ( is_on_positive_halfspace(l[1], p) ||
	 is_on_positive_halfspace(l[1], r) ) {
      num_oriented++;
    }

    if ( is_on_positive_halfspace(l[2], p) ||
	 is_on_positive_halfspace(l[2], q) ) {
      num_oriented++;
    }

    return ( num_oriented >= 2 );
  }

  void
  orient_lines(const Segment& p, const Segment& q, const Segment& r,
		RT a[], RT b[], RT c[]) const 
  {
    compute_supporting_line(p, a[0], b[0], c[0]);
    compute_supporting_line(q, a[1], b[1], c[1]);
    compute_supporting_line(r, a[2], b[2], c[2]);
    
    Line l[3] = {Line(a[0], b[0], c[0]), Line(a[1], b[1], c[1]), 
		 Line(a[2], b[2], c[2])};

    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(l[0], q) ||
	 is_on_positive_halfspace(l[0], r) ) {
      is_oriented[0] = true;
    } else {
      l[0] = l[0].opposite();
      if ( is_on_positive_halfspace(l[0], q) ||
	   is_on_positive_halfspace(l[0], r) ) {
	is_oriented[0] = true;
      } else {
	l[0] = l[0].opposite();
      }
    }

    if ( is_on_positive_halfspace(l[1], p) ||
	 is_on_positive_halfspace(l[1], r) ) {
      is_oriented[1] = true;
    } else {
      l[1] = l[1].opposite();
      if ( is_on_positive_halfspace(l[1], p) ||
	   is_on_positive_halfspace(l[1], r) ) {
	is_oriented[1] = true;
      } else {
	l[1] = l[1].opposite();
      }
    }

    if ( is_on_positive_halfspace(l[2], p) ||
	 is_on_positive_halfspace(l[2], q) ) {
      is_oriented[2] = true;
    } else {
      l[2] = l[2].opposite();
      if ( is_on_positive_halfspace(l[2], p) ||
	   is_on_positive_halfspace(l[2], q) ) {
	is_oriented[2] = true;
      } else {
	l[2] = l[2].opposite();
      }
    }

    if ( is_oriented[0] && is_oriented[1] && is_oriented[2] ) {
      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    int i_no(-1);
    for (int i = 0; i < 3; i++) {
      if ( !is_oriented[i] ) {
	i_no = i;
	CGAL_assertion( is_oriented[(i+1)%3] && is_oriented[(i+2)%3] );
	break;
      }
    }

    CGAL_assertion( i_no != -1 );

    RT  d[3];
    for (int i = 0; i < 3; i++) {
      d[i] = CGAL_NTS square(a[i]) + CGAL_NTS square(b[i]);
    }

    RT z[3];
    for (int i = 0; i < 3; i++) {
      z[i] = l[(i+1)%3].a() * l[(i+2)%3].b()
	- l[(i+2)%3].a() * l[(i+1)%3].b();
    }

    
#if 1
    Sqrt_1 Zero(RT(0), RT(0), d[0]);
    Sqrt_1 sqrt_D0(RT(0), RT(1), d[0]);

    Sqrt_1 D1 = d[1] + Zero;
    Sqrt_1 D2 = d[2] + Zero;

    Sqrt_3 vz(z[0] * sqrt_D0, z[1] + Zero, z[2] + Zero, Zero, D1, D2);

    Sign s_minus_vz = CGAL_NTS sign(vz);
#else
#if 0
    Sign s_minus_vz =
      sign_a_x_sqrt_d_plus_b_x_sqrt_e_plus_c_x_sqrt_f(z[0], z[1], z[2],
						      d[0], d[1], d[2]);
#else
    Sign s_minus_vz = CGAL_NTS sign( z[0] * CGAL_NTS sqrt(d[0]) +
				     z[1] * CGAL_NTS sqrt(d[1]) +
				     z[2] * CGAL_NTS sqrt(d[2]) );
#endif
#endif
    CGAL_assertion( s_minus_vz != ZERO );

    if ( s_minus_vz == NEGATIVE ) {
      l[i_no] = l[i_no].opposite();

      // all these from here...
      z[(i_no+1)%3] = -z[(i_no+1)%3];
      z[(i_no+2)%3] = -z[(i_no+2)%3];


#if 1
      vz =
	Sqrt_3(z[0] * sqrt_D0, z[1] + Zero, z[2] + Zero, Zero, D1, D2);

      s_minus_vz = CGAL_NTS sign(vz);
#else
#if 0
      s_minus_vz =
	sign_a_x_sqrt_d_plus_b_x_sqrt_e_plus_c_x_sqrt_f(z[0], z[1], z[2],
							d[0], d[1], d[2]);
#else
      s_minus_vz = CGAL_NTS sign( z[0] * CGAL_NTS sqrt(d[0]) +
				  z[1] * CGAL_NTS sqrt(d[1]) +
				  z[2] * CGAL_NTS sqrt(d[2]) );
#endif
#endif
      CGAL_assertion( s_minus_vz == POSITIVE );
      // ... up to here are optional; they are computed just for
      // verification.

      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    // now we have to check if the other orientation of l[i_no]
    // corresponds to a CCW triangle as well.
    z[(i_no+1)%3] = -z[(i_no+1)%3];
    z[(i_no+2)%3] = -z[(i_no+2)%3];

#if 1
    vz = Sqrt_3(z[0] * sqrt_D0, z[1] + Zero, z[2] + Zero, Zero, D1, D2);

    Sign s_minus_vz_2 = CGAL_NTS sign(vz);
#else
#if 0
    Sign s_minus_vz_2 =
      sign_a_x_sqrt_d_plus_b_x_sqrt_e_plus_c_x_sqrt_f(z[0], z[1], z[2],
						      d[0], d[1], d[2]);
#else
    Sign s_minus_vz_2 = CGAL_NTS sign( z[0] * CGAL_NTS sqrt(d[0]) +
				       z[1] * CGAL_NTS sqrt(d[1]) +
				       z[2] * CGAL_NTS sqrt(d[2]) );
#endif
#endif
    CGAL_assertion( s_minus_vz_2 != ZERO );

    if ( s_minus_vz_2 == NEGATIVE ) {
      // the other orientation does not correspond to a CCW triangle.
      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    // now compute the Voronoi vertex;
    for (int i = 0; i < 3; i++) {
      a[i] = l[i].a();
      b[i] = l[i].b();
      c[i] = l[i].c();
    }

    RT x[3], y[3], w[3];
    for (int i = 0; i < 3; i++) {
      x[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      y[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      w[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
    }

#if 0
    RT vx = RT(0), vy = RT(0), vw = RT(0);

    for (int i = 0; i < 3; i++) {
      RT D = CGAL_NTS sqrt(d[i]);
      vx += x[i] * D;
      vy += y[i] * D;
      vw += w[i] * D;
    }

    RT dist = vw * ( a[(i_no+1)%3] * vx + b[(i_no+1)%3] * vy +
		     c[(i_no+1)%3] * vw );
    
#else
    Sqrt_3 dist;
    {
      Sqrt_1 Zero(RT(0), RT(0), d[0]);
      Sqrt_1 sqrt_D0(RT(0), RT(1), d[0]);

      Sqrt_1 D1 = d[1] + Zero;
      Sqrt_1 D2 = d[2] + Zero;

      Sqrt_3 vx, vy, vw;

      vx = Sqrt_3(x[0] * sqrt_D0, x[1], x[2], Zero, D1, D2);
      vy = Sqrt_3(y[0] * sqrt_D0, y[1], y[2], Zero, D1, D2);
      vw = Sqrt_3(w[0] * sqrt_D0, w[1], w[2], Zero, D1, D2);

      Sqrt_1 a1(a[(i_no+1)%3], RT(0), d[0]);
      Sqrt_1 b1(b[(i_no+1)%3], RT(0), d[0]);
      Sqrt_1 c1(c[(i_no+1)%3], RT(0), d[0]);

      dist = vw * ( a1 * vx + b1 * vy +	c1 * vw );
    }
#endif


    Sign sgn_dist = CGAL_NTS sign(dist);

    CGAL_assertion( sgn_dist != ZERO );

    if ( sgn_dist == NEGATIVE ) {
      a[i_no] = -a[i_no];
      b[i_no] = -b[i_no];
      c[i_no] = -c[i_no];
    }
  }



#if 0
  void
  orient_lines1(const Segment& p, const Segment& q, const Segment& r,
		RT a[], RT b[], RT c[]) const
  {
    compute_supporting_line(p, a[0], b[0], c[0]);
    compute_supporting_line(q, a[1], b[1], c[1]);
    compute_supporting_line(r, a[2], b[2], c[2]);
    
    Line l[3] = {Line(a[0], b[0], c[0]), Line(a[1], b[1], c[1]), 
		 Line(a[2], b[2], c[2])};

    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(l[0], q) ||
	 is_on_positive_halfspace(l[0], r) ) {
      is_oriented[0] = true;
    } else {
      l[0] = l[0].opposite();
      if ( is_on_positive_halfspace(l[0], q) ||
	   is_on_positive_halfspace(l[0], r) ) {
	is_oriented[0] = true;
      } else {
	l[0] = l[0].opposite();
      }
    }

    if ( is_on_positive_halfspace(l[1], p) ||
	 is_on_positive_halfspace(l[1], r) ) {
      is_oriented[1] = true;
    } else {
      l[1] = l[1].opposite();
      if ( is_on_positive_halfspace(l[1], p) ||
	   is_on_positive_halfspace(l[1], r) ) {
	is_oriented[1] = true;
      } else {
	l[1] = l[1].opposite();
      }
    }

    if ( is_on_positive_halfspace(l[2], p) ||
	 is_on_positive_halfspace(l[2], q) ) {
      is_oriented[2] = true;
    } else {
      l[2] = l[2].opposite();
      if ( is_on_positive_halfspace(l[2], p) ||
	   is_on_positive_halfspace(l[2], q) ) {
	is_oriented[2] = true;
      } else {
	l[2] = l[2].opposite();
      }
    }

    if ( is_oriented[0] && is_oriented[1] && is_oriented[2] ) {
      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    int i_no(-1);
    for (int i = 0; i < 3; i++) {
      if ( !is_oriented[i] ) {
	i_no = i;
	CGAL_assertion( is_oriented[(i+1)%3] && is_oriented[(i+2)%3] );
	break;
      }
    }

    CGAL_assertion( i_no != -1 );

    RT  d[3];
    for (int i = 0; i < 3; i++) {
      d[i] = CGAL_NTS square(a[i]) + CGAL_NTS square(b[i]);
    }

    RT z[3];
    for (int i = 0; i < 3; i++) {
      z[i] = l[(i+1)%3].a() * l[(i+2)%3].b()
	- l[(i+2)%3].a() * l[(i+1)%3].b();
    }

    Sign s_minus_vz = CGAL_NTS sign( z[0] * CGAL_NTS sqrt(d[0]) +
				     z[1] * CGAL_NTS sqrt(d[1]) +
				     z[2] * CGAL_NTS sqrt(d[2]) );

    CGAL_assertion( s_minus_vz != ZERO );

    if ( s_minus_vz == NEGATIVE ) {
      l[i_no] = l[i_no].opposite();
    }
    for (int i = 0; i < 3; i++) {
      a[i] = l[i].a();
      b[i] = l[i].b();
      c[i] = l[i].c();
    }
  }

  void
  orient_lines2(const Segment& p, const Segment& q, const Segment& r,
		RT a[], RT b[], RT c[]) const
  {
    compute_supporting_line(p, a[0], b[0], c[0]);
    compute_supporting_line(q, a[1], b[1], c[1]);
    compute_supporting_line(r, a[2], b[2], c[2]);

    RT  d[3];
    for (int i = 0; i < 3; i++) {
      d[i] = CGAL_NTS square(a[i]) + CGAL_NTS square(b[i]);
    }

    int s1[] = {1,  1,  1,  1, -1, -1, -1, -1};
    int s2[] = {1,  1, -1, -1,  1,  1, -1, -1};
    int s3[] = {1, -1,  1, -1,  1, -1,  1, -1};
    bool is_pos[] = {false,false,false,false,false,false,false,false};
    bool is_ok[] = {false,false,false,false,false,false,false,false};
    int num_pos(0);

    int num_cases = 8;

    for (int k = 0; k < num_cases; k++) {
      int sgn[3] = {s1[k], s2[k], s3[k]};
      RT al[3], bl[3], cl[3];
      RT z[3];

      for (int i = 0; i < 3; i++) {
	if ( sgn[i] > 0 ) {
	  al[i] = a[i];
	  bl[i] = b[i];
	  cl[i] = c[i];
	} else {
	  al[i] = -a[i];
	  bl[i] = -b[i];
	  cl[i] = -c[i];
	}
      }

      for (int i = 0; i < 3; i++) {
	z[i] = al[(i+1)%3] * bl[(i+2)%3] - al[(i+2)%3] * bl[(i+1)%3];
      }

#if 0
      Sign s_minus_vz =
	sign_a_x_sqrt_d_plus_b_x_sqrt_e_plus_c_x_sqrt_f(z[0], z[1], z[2],
							d[0], d[1], d[2]);
#else
      Sign s_minus_vz =
	CGAL_NTS sign( z[0] * CGAL_NTS sqrt(d[0]) +
		       z[1] * CGAL_NTS sqrt(d[1]) +
		       z[2] * CGAL_NTS sqrt(d[2]) );
#endif
      is_pos[k] = ( s_minus_vz == POSITIVE );
      if ( s_minus_vz == POSITIVE ) {
	num_pos++;
	is_ok[k] = is_consistent(p, q, r, al, bl, cl);
      }
    }

    if ( num_pos == 1 ) {
      for (int k = 0; k < num_cases; k++) {
	if ( is_pos[k] ) {
	  is_ok[k] = true;
	  for (int i = 1; i < num_cases; i++) {
	    is_ok[(k+i)%num_cases] = false;
	  }
	  break;
	}
      }
    } else {
      //      CGAL_assertion( num_pos == 3 );
      for (int k = 0; k < num_cases; k++) {
	if ( is_pos[k] && is_ok[k] ) {
	  for (int i = 1; i < num_cases; i++) {
	    //	    CGAL_assertion( !is_pos[(k+i)%num_cases] |
	    //	    !is_ok[(k+i)%num_cases] );
	  }
	}
      }
    }

    for (int k = 0; k < num_cases; k++) {
      if ( is_pos[k] && is_ok[k] ) {
	int sgn[3] = {s1[k], s2[k], s3[k]};
	for (int i = 0; i < 3; i++) {
	  if ( sgn[i] < 0 ) {
	    a[i] = -a[i];
	    b[i] = -b[i];
	    c[i] = -c[i];
	  }
	}
	return;
      } // end of outer-if
    }

  }
#endif


  void
  compute_sss(const Segment& p, const Segment& q, const Segment& r)
  {
    v_type = SSS;

    RT a[3], b[3], c[3];

    orient_lines(p, q, r, a, b, c);

    for (int i = 0; i < 3; i++) {
      cx[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      cy[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      cz[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
      D[i] = CGAL_NTS square(a[i]) + CGAL_NTS square(b[i]);
    }
  }

  //--------------------------------------------------------------------------

  void
  compute_vertex(const Site& s1, const Site& s2, const Site& s3)
  {
    if ( s1.is_point() && s2.is_point() && s3.is_point() ) {
      compute_ppp(s1.point(), s2.point(), s3.point());

    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      compute_vertex(s2, s3, s1);
      pps_idx = 1;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
      pps_idx = 2;

    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      
      compute_pps(s1.point(), s2.point(), s3.segment());
      pps_idx = 0;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      compute_pss(s1.point(), s2.segment(), s3.segment());
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      compute_vertex(s2, s3, s1);
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
    } else {
      compute_sss(s1.segment(), s2.segment(), s3.segment());
    }

  }


  //--------------------------------------------------------------------------

  bool are_identical(const Point& p, const Point& q) const
  {
    return (p == q);
  }

  //--------------------------------------------------------------------------

  bool is_endpoint_of(const Point& p, const Segment& s) const
  {
    return ( are_identical(p, s.source()) ||
	     are_identical(p, s.target()) );
  }
  

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           the orientation test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  template<class Type>
  Orientation
  orientation(const Line& l, Type, Sqrt_field_tag) const
  {
    Sign s = CGAL_NTS sign(l.a() * x() + l.b() * y() + l.c());

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }
  
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line& l, PPP_Type, Ring_tag) const
  {
    Sign s_uz = CGAL_NTS sign(uz);
    Sign s_l = CGAL_NTS sign(l.a() * ux + l.b() * uy + l.c() * uz);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }

  
  Orientation
  orientation(const Line& l, PPS_Type, Ring_tag) const
  {
    Sqrt_1 vx(J + p_ref().x() * X, RT(-2) * y_, S);
    Sqrt_1 vy(I + p_ref().y() * X, RT( 2) * x_, S);
    Sqrt_1 vz(X, RT(0), S);

    Sign s_uz = CGAL_NTS sign(X);
    Sign s_l = CGAL_NTS sign(l.a() * vx + l.b() * vy + l.c() * vz);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }

  Orientation
  orientation(const Line& l, PSS_Type, Ring_tag) const
  {
    RT B = a1a2 + b1b2;
    Sqrt_1 vz1(-B, RT(1), D1D2);

    RT A = a1a2 - b1b2;
    Sqrt_1 u1( c1c2 * A, c1c2, D1D2);
    Sqrt_1 u2(-c1c2 * A, c1c2, D1D2);


    int is(_s);
    int ir(_r);
    RT ss(is);
    RT rr(ir);

    Sqrt_1 Zero(RT(0), RT(0), D1D2);
    Sqrt_1 sigma(ss, RT(0), D1D2);
    Sqrt_1 rho(rr, RT(0), D1D2);

    Sqrt_1 J1(J, RT(0), D1D2);
    Sqrt_1 I1(I, RT(0), D1D2);

    Sqrt_3 vx(J1 + p_ref().x() * vz1, sigma, Zero, Zero, u1, u2);
    Sqrt_3 vy(I1 + p_ref().y() * vz1, Zero, -rho, Zero, u1, u2);

    Sqrt_3 vz(vz1, Zero, Zero, Zero, u1, u2);

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;
    Sqrt_1 c = l.c() + Zero;

    Sign s_uz = CGAL_NTS sign(vz1);
    Sign s_l = CGAL_NTS sign(a * vx + b * vy + c * vz);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }
    
  Orientation
  orientation(const Line& l, SSS_Type, Ring_tag) const
  {
    Sqrt_1 Zero(RT(0), RT(0), D[0]);
    Sqrt_1 sqrt_D0(RT(0), RT(1), D[0]);

    Sqrt_1 D1 = D[1] + Zero;
    Sqrt_1 D2 = D[2] + Zero;

    Sqrt_3 vx(cx[0] * sqrt_D0, cx[1] + Zero, cx[2] + Zero, Zero, D1, D2);
    Sqrt_3 vy(cy[0] * sqrt_D0, cy[1] + Zero, cy[2] + Zero, Zero, D1, D2);
    Sqrt_3 vz(cz[0] * sqrt_D0, cz[1] + Zero, cz[2] + Zero, Zero, D1, D2);


    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;
    Sqrt_1 c = l.c() + Zero;

    Sign s_uz = CGAL_NTS sign(vz);
    Sign s_l = CGAL_NTS sign(a * vx + b * vy + c * vz);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                              the incircle test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a point
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Sign
  check_easy_degeneracies(const Point& t, PPS_Type,
			  bool& use_result) const
  {
    use_result = false;
    if (  ( p_.is_point() && are_identical(p_.point(),t) ) ||
	  ( q_.is_point() && are_identical(q_.point(),t) ) ||
	  ( r_.is_point() && are_identical(r_.point(),t) )  ) {
      use_result = true;
      return ZERO;
    }

    if (  ( p_.is_segment() && is_endpoint_of(t, p_.segment()) ) || 
	  ( q_.is_segment() && is_endpoint_of(t, q_.segment()) ) ||
	  ( r_.is_segment() && is_endpoint_of(t, r_.segment()) )  ) {
      use_result = true;
      return POSITIVE;
    }

    return ZERO;
  }

  Sign
  check_easy_degeneracies(const Point& t, PSS_Type,
			  bool& use_result) const
  {
    return check_easy_degeneracies(t, PPS_Type(), use_result);
  }

  Sign
  check_easy_degeneracies(const Point& t, SSS_Type,
			  bool& use_result) const
  {
    use_result = false;

    // ADD THE CASES WHERE t IS AN ENDPOINT OF ONE OF THE SEGMENTS
    return ZERO;
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Point& t, PPP_Type) const
  {
    Oriented_side os =
      side_of_oriented_circle(p_.point(), q_.point(), r_.point(), t);
    if ( os == ON_POSITIVE_SIDE ) { return NEGATIVE; }
    if ( os == ON_NEGATIVE_SIDE ) { return POSITIVE; }
    return ZERO;
  }

  //--------------------------------------------------------------------------

  template<class Type>
  Sign incircle(const Point& t, Type type, Sqrt_field_tag) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    FT r2 = squared_radius();

    FT d2 = CGAL_NTS square(x() - t.x()) +
      CGAL_NTS square(y() - t.y());

    return Sign( CGAL_NTS compare(d2, r2) );
  }

  //--------------------------------------------------------------------------

  
  Sign incircle(const Point& t, PPS_Type type, Ring_tag) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    Sqrt_1 vx(J, RT(-2) * y_, S);
    Sqrt_1 vy(I, RT( 2) * x_, S);

    Sqrt_1 Rs = CGAL_NTS square(vx) + CGAL_NTS square(vy);

    RT dx = p_ref().x() - t.x();
    RT dy = p_ref().y() - t.y();

    Sqrt_1 Rs1 =
      CGAL_NTS square(vx + RT(dx * X)) + CGAL_NTS square(vy + RT(dy * X));

    return CGAL_NTS sign(Rs1 - Rs);
  }

  Sign incircle(const Point& t, PSS_Type type, Ring_tag) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    RT B = a1a2 + b1b2;
    Sqrt_1 vz(-B, RT(1), D1D2);

    RT A = a1a2 - b1b2;
    Sqrt_1 u1( c1c2 * A, c1c2, D1D2);
    Sqrt_1 u2(-c1c2 * A, c1c2, D1D2);

#if 1
    Sqrt_3 vx, vy;

    Sqrt_1 Zero, sigma, rho;

    Zero = Sqrt_1(RT(0), RT(0), D1D2);
    sigma = Sqrt_1(RT(int(_s)), RT(0), D1D2);
    rho = Sqrt_1(RT(int(_r)), RT(0), D1D2);

    vx = Sqrt_3(Sqrt_1(J, RT(0), D1D2), sigma, Zero, Zero, u1, u2);
    vy = Sqrt_3(Sqrt_1(I, RT(0), D1D2), Zero, -rho, Zero, u1, u2);
#else
    Sqrt_3 vx(Sqrt_1(J), Sqrt_1( int(_s)), Sqrt_1(0), Sqrt_1(0), u1, u2);
    Sqrt_3 vy(Sqrt_1(I), Sqrt_1(0), Sqrt_1(-int(_r)), Sqrt_1(0), u1, u2);
#endif

    Sqrt_3 Rs = CGAL_NTS square(vx) + CGAL_NTS square(vy);


#if 0
    Sqrt_3 dx((p_ref().x() - t.x()) * vz, Zero, Zero, Zero, u1, u2);
    Sqrt_3 dy((p_ref().y() - t.y()) * vz, Zero, Zero, Zero, u1, u2);
#else
    Sqrt_1 dxvz = RT((p_ref().x() - t.x())) * vz;
    Sqrt_1 dyvz = RT((p_ref().y() - t.y())) * vz;
    Sqrt_3 dx(dxvz, Zero, Zero, Zero, u1, u2);
    Sqrt_3 dy(dyvz, Zero, Zero, Zero, u1, u2);
#endif

#if 0
    Sqrt_3 vx_dx = vx + dx;
    Sqrt_3 vy_dy = vy + dy;
    Sqrt_3 Rs1 = CGAL_NTS square(vx_dx) + CGAL_NTS square(vy_dy);
#else
    Sqrt_3 Rs1 = CGAL_NTS square(vx + dx) + CGAL_NTS square(vy + dy);
#endif

    //    Sign s_vz = CGAL_NTS sign(vz);
    Sign s_Q = CGAL_NTS sign(Rs1 - Rs);

    //    return Sign(s_vz * s_Q);
    return s_Q;
  }


  Sign incircle(const Point& t, SSS_Type, Ring_tag) const
  {
    Sqrt_1 Zero(RT(0), RT(0), D[0]);
    Sqrt_1 One = Zero + RT(1);
    Sqrt_1 sqrt_D0(RT(0), RT(1), D[0]);

    Sqrt_1 D1 = D[1] + Zero;
    Sqrt_1 D2 = D[2] + Zero;

    Sqrt_3 vx(cx[0] * sqrt_D0, cx[1] + Zero, cx[2] + Zero, Zero, D1, D2);
    Sqrt_3 vy(cy[0] * sqrt_D0, cy[1] + Zero, cy[2] + Zero, Zero, D1, D2);
    Sqrt_3 vz(cz[0] * sqrt_D0, cz[1] + Zero, cz[2] + Zero, Zero, D1, D2);

    RT a1, b1, c1;
    compute_supporting_line(p_.segment(), a1, b1, c1);

    Sqrt_1 a = a1 + Zero;
    Sqrt_1 b = b1 + Zero;
    Sqrt_1 c = c1 + Zero;

    Sqrt_1 Ns = CGAL_NTS square(a) + CGAL_NTS square(b);
    Sqrt_3 Ls = CGAL_NTS square(a * vx + b * vy + c * vz);

    Sqrt_1 tx = t.x() + Zero;
    Sqrt_1 ty = t.y() + Zero;

    Sqrt_3 R1s = CGAL_NTS square(vx - tx * vz)
      + CGAL_NTS square(vy - ty * vz);

    return CGAL_NTS sign(R1s * Ns - Ls);
  }



  //--------------------------------------------------------------------------

  template<class Method_tag>
  Sign incircle(const Point& t, Method_tag tag) const 
  {
    if ( is_degenerate_Voronoi_circle() ) {
      return POSITIVE;
    }

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle(t, PPP_Type());
      break;
    case PPS:
      s = incircle(t, PPS_Type(), tag);
      break;
    case PSS:
      s = incircle(t, PSS_Type(), tag);
      break;
    case SSS:
      s = incircle(t, SSS_Type(), tag);
      break;
    }

    return s;
  }

  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a segment
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Oriented_side
  oriented_side(const Line& l, const Point& p, PPP_Type, Ring_tag) const
  {
    Sign s_uz = CGAL_NTS sign(uz);

    RT px = uz * p.x() - ux;
    RT py = uz * p.y() - uy;

    Sign s1 = CGAL_NTS sign(l.b() * px - l.a() * py);

    Sign s = Sign(s_uz * s1);

    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }

  Oriented_side
  oriented_side(const Line& l, const Point& p, PPS_Type, Ring_tag) const
  {
    Sqrt_1 vx1(J + (p_ref().x() - p.x()) * X, RT(-2) * y_, S);
    Sqrt_1 vy1(I + (p_ref().y() - p.y()) * X, RT( 2) * x_, S);

    Sign s = Sign(CGAL_NTS sign(X) *
		  CGAL_NTS sign(vy1 * l.a() - vx1 * l.b())
		  );

    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }


  Oriented_side
  oriented_side(const Line& l, const Point& p, PSS_Type, Ring_tag) const
  {
    RT B = a1a2 + b1b2;
    Sqrt_1 vz(-B, RT(1), D1D2);

    RT A = a1a2 - b1b2;
    Sqrt_1 u1( c1c2 * A, c1c2, D1D2);
    Sqrt_1 u2(-c1c2 * A, c1c2, D1D2);

#if 0
    Sqrt_3 vx, vy;

    Sqrt_1 Zero, sigma, rho;

    Zero = Sqrt_1(RT(0), RT(0), D1D2);
    sigma = Sqrt_1(RT(int(_s)), RT(0), D1D2);
    rho = Sqrt_1(RT(int(_r)), RT(0), D1D2);

    vx = Sqrt_3(Sqrt_1(J, RT(0), D1D2), sigma, Zero, Zero, u1, u2);
    vy = Sqrt_3(Sqrt_1(I, RT(0), D1D2), Zero, -rho, Zero, u1, u2);
#else

    int is(_s);
    int ir(_r);
    RT ss(is);
    RT rr(ir);

    Sqrt_1 Zero(RT(0), RT(0), D1D2);
    Sqrt_1 sigma(ss, RT(0), D1D2);
    Sqrt_1 rho(rr, RT(0), D1D2);

    Sqrt_3 vx(Sqrt_1(J, RT(0), D1D2), sigma, Zero, Zero, u1, u2);
    Sqrt_3 vy(Sqrt_1(I, RT(0), D1D2), Zero, -rho, Zero, u1, u2);

#endif

    RT dxx = p_ref().x() - p.x();
    RT dyy = p_ref().y() - p.y();

    Sqrt_3 dx(dxx * vz, Zero, Zero, Zero, u1, u2);
    Sqrt_3 dy(dyy * vz, Zero, Zero, Zero, u1, u2);

    Sqrt_3 vx1 = vx + dx;
    Sqrt_3 vy1 = vy + dy;

    Sqrt_1 a(l.a(), RT(0), D1D2);
    Sqrt_1 b(l.b(), RT(0), D1D2);

    Sign s = Sign(CGAL_NTS sign(vz) *
		  CGAL_NTS sign(a * vy1 - b * vx1)
		  );


    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }

  Oriented_side
  oriented_side(const Line& l, const Point& p, SSS_Type, Ring_tag) const
  {
    Sqrt_1 Zero(RT(0), RT(0), D[0]);
    //    Sqrt_1 One = Zero + RT(1);
    Sqrt_1 sqrt_D0(RT(0), RT(1), D[0]);

    Sqrt_1 D1 = D[1] + Zero;
    Sqrt_1 D2 = D[2] + Zero;

    Sqrt_3 vx(cx[0] * sqrt_D0, cx[1] + Zero, cx[2] + Zero, Zero, D1, D2);
    Sqrt_3 vy(cy[0] * sqrt_D0, cy[1] + Zero, cy[2] + Zero, Zero, D1, D2);
    Sqrt_3 vz(cz[0] * sqrt_D0, cz[1] + Zero, cz[2] + Zero, Zero, D1, D2);

    Sqrt_1 px = p.x() + Zero;
    Sqrt_1 py = p.y() + Zero;

    Sqrt_3 vx1 = vx - px * vz;
    Sqrt_3 vy1 = vy - py * vz;

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;

    Sign s = Sign(CGAL_NTS sign(vz) *
		  CGAL_NTS sign(a * vy1 - b * vx1)
		  );

    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }
  

  //--------------------------------------------------------------------------

  
  Sign incircle(const Line& l, PPP_Type, Ring_tag) const
  {
    RT a1 = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());
    RT a2 = CGAL_NTS square(ux - p_ref().x() * uz) +
      CGAL_NTS square(uy - p_ref().y() * uz);

    RT a3 = CGAL_NTS square(l.a() * ux + l.b() * uy + l.c() * uz);

    Comparison_result cr = CGAL_NTS compare(a3, a1 * a2);

    if ( cr == LARGER ) { return POSITIVE; }
    if ( cr == SMALLER ) { return NEGATIVE; }
    return ZERO;
  }

  Sign incircle(const Line& l, PPS_Type, Ring_tag) const
  {
    Sqrt_1 vx(J, RT(-2) * y_, S);
    Sqrt_1 vy(I, RT( 2) * x_, S);

    Sqrt_1 Rs = CGAL_NTS square(vx) + CGAL_NTS square(vy);

    Sqrt_1 vx1 = vx + RT(p_ref().x() * X);
    Sqrt_1 vy1 = vy + RT(p_ref().y() * X);
    
    RT Ns = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

    Sqrt_1 Ls = CGAL_NTS square(RT(l.a()) * vx1
				+ RT(l.b()) * vy1
				+ RT(l.c() * X));

    return CGAL_NTS sign(Ls - Rs * Ns);
  }


  Sign incircle(const Line& l, PSS_Type, Ring_tag) const
  {
    RT B = a1a2 + b1b2;
    Sqrt_1 vz(-B, RT(1), D1D2);

    RT A = a1a2 - b1b2;
    Sqrt_1 u1( c1c2 * A, c1c2, D1D2);
    Sqrt_1 u2(-c1c2 * A, c1c2, D1D2);

#if 1
    Sqrt_3 vx, vy;

    Sqrt_1 Zero, sigma, rho;

    Zero = Sqrt_1(RT(0), RT(0), D1D2);
    sigma = Sqrt_1(RT(int(_s)), RT(0), D1D2);
    rho = Sqrt_1(RT(int(_r)), RT(0), D1D2);

    vx = Sqrt_3(Sqrt_1(J, RT(0), D1D2), sigma, Zero, Zero, u1, u2);
    vy = Sqrt_3(Sqrt_1(I, RT(0), D1D2), Zero, -rho, Zero, u1, u2);
#else
    Sqrt_3 vx(Sqrt_1(J), Sqrt_1( int(_s)), Sqrt_1(0), Sqrt_1(0), u1, u2);
    Sqrt_3 vy(Sqrt_1(I), Sqrt_1(0), Sqrt_1(-int(_r)), Sqrt_1(0), u1, u2);
#endif

    Sqrt_3 Rs = CGAL_NTS square(vx) + CGAL_NTS square(vy);

    Sqrt_3 dx(p_ref().x() * vz, Zero, Zero, Zero, u1, u2);
    Sqrt_3 dy(p_ref().y() * vz, Zero, Zero, Zero, u1, u2);

    Sqrt_3 vx1 = vx + dx;
    Sqrt_3 vy1 = vy + dy;
    Sqrt_3 vz1(vz, Zero, Zero, Zero, u1, u2);

    Sqrt_1 a(l.a(), RT(0), D1D2);
    Sqrt_1 b(l.b(), RT(0), D1D2);
    Sqrt_1 c(l.c(), RT(0), D1D2);

    Sqrt_1 Ns = CGAL_NTS square(a) + CGAL_NTS square(b);

    Sqrt_3 Ls = CGAL_NTS square(a * vx1 + b * vy1 + c * vz1);

    return CGAL_NTS sign(Ls - Rs * Ns);
  }


  Sign incircle(const Line& l, SSS_Type, Ring_tag) const
  {
    Sqrt_1 Zero(RT(0), RT(0), D[0]);
    Sqrt_1 One = Zero + RT(1);
    Sqrt_1 sqrt_D0(RT(0), RT(1), D[0]);

    Sqrt_1 D1 = D[1] + Zero;
    Sqrt_1 D2 = D[2] + Zero;

    Sqrt_3 vx(cx[0] * sqrt_D0, cx[1] + Zero, cx[2] + Zero, Zero, D1, D2);
    Sqrt_3 vy(cy[0] * sqrt_D0, cy[1] + Zero, cy[2] + Zero, Zero, D1, D2);
    Sqrt_3 vz(cz[0] * sqrt_D0, cz[1] + Zero, cz[2] + Zero, Zero, D1, D2);

    RT a1, b1, c1;
    compute_supporting_line(p_.segment(), a1, b1, c1);

    Sqrt_1 a = a1 + Zero;
    Sqrt_1 b = b1 + Zero;
    Sqrt_1 c = c1 + Zero;

    Sqrt_1 Ns = CGAL_NTS square(a) + CGAL_NTS square(b);

    Sqrt_1 la = l.a() + Zero;
    Sqrt_1 lb = l.b() + Zero;
    Sqrt_1 lc = l.c() + Zero;

    Sqrt_1 Ns1 = CGAL_NTS square(la) + CGAL_NTS square(lb);

    Sqrt_3 Ls = CGAL_NTS square(a * vx + b * vy + c * vz);

    Sqrt_3 Ls1 =
      CGAL_NTS square(la * vx + lb * vy + lc * vz);

    return CGAL_NTS sign(Ls1 * Ns - Ls * Ns1);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  template<class Type>
  Oriented_side
  oriented_side(const Line& l, const Point& p, Type, Sqrt_field_tag) const
  {
    Line l1(l.b(), -l.a(), l.a() * y() - l.b() * x());

    return l1.oriented_side(p);
  }


  template<class Type>
  Sign incircle(const Line& l, Type, Sqrt_field_tag) const
  {
    RT r2 = squared_radius();

    RT n2 = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

    RT d2 = CGAL_NTS square(l.a() * x() + l.b() * y() + l.c());
    //    RT d2 = CGAL_NTS square(l.a() * x() + l.b() * y() + l.c()) / n2;

    return Sign( CGAL_NTS compare(d2, r2 * n2) );
    //    return Sign( CGAL_NTS compare(d2, r2) );
  }



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  template<class Type, class Method_tag>
  Sign incircle(const Segment& t, Type type, Method_tag tag) const
  {
    if ( v_type == PPP || v_type == PPS ) {
      if (  p_.is_point() && q_.is_point() &&
	    is_endpoint_of(p_.point(), t) &&
	    is_endpoint_of(q_.point(), t)  ) {
	return NEGATIVE;
      }

      if (  p_.is_point() && r_.is_point() &&
	    is_endpoint_of(p_.point(), t) &&
	    is_endpoint_of(r_.point(), t)  ){
	return NEGATIVE;
      }

      if (  q_.is_point() && r_.is_point() &&
	    is_endpoint_of(q_.point(), t) &&
	    is_endpoint_of(r_.point(), t)  ){
	return NEGATIVE;
      }
    }

    Sign d1, d2;
    if (  ( p_.is_point() && are_identical(p_.point(), t.source()) ) ||
	  ( q_.is_point() && are_identical(q_.point(), t.source()) ) ||
	  ( r_.is_point() && are_identical(r_.point(), t.source()) )  ) {
      d1 = ZERO;
    } else {
      d1 = incircle(t.source(), tag);
      //      d1 = incircle(t.source(), type, tag);
    }
    if ( d1 == NEGATIVE ) { return NEGATIVE; }

    if (  ( p_.is_point() && are_identical(p_.point(), t.target()) ) ||
	  ( q_.is_point() && are_identical(q_.point(), t.target()) ) ||
	  ( r_.is_point() && are_identical(r_.point(), t.target()) )  ) {
      d2 = ZERO;
    } else {
      d2 = incircle(t.target(), tag);
      //      d2 = incircle(t.target(), type, tag);
    }
    if ( d2 == NEGATIVE ) { return NEGATIVE; }

    RT a, b, c;
    compute_supporting_line(t, a, b, c);

    // MK:: THIS IS REALLY IMPORTANT BECAUSE IT DEPENDS ON HOW LINES
    //    ARE IMPLEMENTED *******************************************
    //      change from passing a line to passing a triple of numbers
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Line l(a, b, c);
    Sign sl = incircle(l, type, tag);

    // this is the old code
    //    if ( sl == NEGATIVE ) { return sl; }
    if ( sl == POSITIVE ) { return sl; }

    Oriented_side os1 = oriented_side(l, t.source(), type, tag);
    Oriented_side os2 = oriented_side(l, t.target(), type, tag);

    if ( sl == ZERO ) {
      if ( (os1 == ON_POSITIVE_SIDE && os2 == ON_NEGATIVE_SIDE) ||
	   (os1 == ON_NEGATIVE_SIDE && os2 == ON_POSITIVE_SIDE) ) {
	return ZERO;
      }
      return POSITIVE;
    }

    return (os1 == os2) ? POSITIVE : NEGATIVE;
  }

  //--------------------------------------------------------------------------

  template<class Method_tag>
  Sign incircle(const Segment& t, Method_tag tag) const 
  {
    if ( is_degenerate_Voronoi_circle() ) {
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle
      if (  !are_identical( p_ref(), t.source() ) &&
	    !are_identical( p_ref(), t.target() )  ) {
	return POSITIVE;
      }

      if (  ( r_.is_point() && are_identical(p_ref(),r_.point()) ) ||
	    ( q_.is_point() && are_identical(p_ref(),q_.point()) ) ||
	    ( p_.is_point() && are_identical(p_ref(),p_.point()) )  ) {
	Point pr;
	Segment sp, sq;
	if ( p_.is_point() ) {
	  CGAL_assertion( q_.is_segment() && r_.is_segment() );
	  pr = p_.point();
	  sp = q_.segment();
	  sq = r_.segment();
	} else if ( q_.is_point() ) {
	  CGAL_assertion( r_.is_segment() && p_.is_segment() );
	  pr = q_.point();
	  sp = r_.segment();
	  sq = p_.segment();
	} else {
	  CGAL_assertion( p_.is_segment() && q_.is_segment() );
	  pr = r_.point();
	  sp = p_.segment();
	  sq = q_.segment();
	}

	Point pq = sq.source(), pp = sp.source(), pt = t.source();

	if ( are_identical(sp.source(),pr) ) { pp = sp.target(); }
	if ( are_identical(sq.source(),pr) ) { pq = sq.target(); }
	if ( are_identical( t.source(),pr) ) { pt =  t.target(); }

	if ( CGAL::orientation(pr, pp, pt) == LEFT_TURN &&
	     CGAL::orientation(pr, pq, pt) == RIGHT_TURN ) {
	  return NEGATIVE;
	}
	return ZERO;
      }
    } // if ( is_degenerate_Voronoi_circle() )

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle(t, PPP_Type(), tag);
      break;
    case PPS:
      s = incircle(t, PPS_Type(), tag);
      break;
    case PSS:
      s = incircle(t, PSS_Type(), tag);
      break;
    case SSS:
      s = incircle(t, SSS_Type(), tag);
      break;
    }

    return s;
  }

  //--------------------------------------------------------------------------
  //  subpredicates for the incircle test
  //--------------------------------------------------------------------------


public:
  bool is_degenerate_Voronoi_circle() const
  {
    if ( v_type != PSS ) { return false; }

    if ( p_.is_point() ) {
      return ( is_endpoint_of(p_.point(), q_.segment()) &&
	       is_endpoint_of(p_.point(), r_.segment()) );
    } else if ( q_.is_point() ) {
      return ( is_endpoint_of(q_.point(), p_.segment()) &&
	       is_endpoint_of(q_.point(), r_.segment()) );
    } else {
      CGAL_assertion( r_.is_point() );
      return ( is_endpoint_of(r_.point(), p_.segment()) &&
	       is_endpoint_of(r_.point(), q_.segment()) );
    }
  }


  //--------------------------------------------------------------------------


public:
  // THIS METHOD SHOULD BE PRIVATE AND SHOULD NOT BE USED IN THE
  // PREDICATES (MAYBE I AM WRONG ON THIS????) IN ANY CASE CALLS TO
  // THIS METHOD SHOULD BE ELIMINATED AS MUCH AS POSSIBLE; ALSO I NEED
  // TO WRITE A SPECIAL VERSION FOR THE RING_TAG CASE
  template<class Method_tag>
  bool is_same_point(const Point& p, Method_tag) const
  {
    Comparison_result res = CGAL_NTS compare(x(), p.x());
    if ( res != EQUAL ) { return false; }

    return (CGAL_NTS compare(y(), p.y()) == EQUAL);
  }

private:

  //--------------------------------------------------------------------------

  bool is_zero_radius(Ring_tag) const
  {
    // to be filled in
    return false;
  }

  bool is_zero_radius(Sqrt_field_tag) const
  {
    return CGAL_NTS is_zero( squared_radius() );
  }

  //--------------------------------------------------------------------------
  //  the reference point (valid if v_type != SSS)
  //--------------------------------------------------------------------------

  const Point& p_ref() const
  {
    CGAL_precondition ( v_type != SSS );

    if ( v_type == PPS ) {
      if ( pps_idx == 0 ) { return p_.point(); }
      if ( pps_idx == 1 ) { return q_.point(); }
      return r_.point();
    }

    if ( p_.is_point() ) {
      return p_.point();
    } else if ( q_.is_point() ) {
      return q_.point();
    } else {
      CGAL_assertion( r_.is_point() );
      return r_.point();
    }
  }


public:
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           access methods
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  inline FT x() const { return hx() / hw(); }
  inline FT y() const { return hy() / hw(); }

  RT hx() const {
    switch (v_type) {
    case PPP:
      {
	return ux;
      }
      break;
    case PPS:
      {
	return J + p_ref().x() * X - FT(2) * y_ * CGAL_NTS sqrt(S);
      }
      break;
    case PSS:
      {
	RT G = c1c2 * (CGAL_NTS sqrt(D1D2) + a1a2 - b1b2);
	return J + p_ref().x() * hw () + RT(int(_s)) * CGAL_NTS sqrt(G);
      }
      break;
    case SSS:
      {
	return ( cx[0] * CGAL_NTS sqrt(D[0]) +
		 cx[1] * CGAL_NTS sqrt(D[1]) +
		 cx[2] * CGAL_NTS sqrt(D[2]) );
      }
      break;
    default:
      return RT(0);
    }
  }

  RT hy() const {
    switch (v_type) {
    case PPP:
      {
	return uy;
      }
      break;
    case PPS:
      {
	return I + p_ref().y() * X + RT(2) * x_ * CGAL_NTS sqrt(S);
      }
      break;
    case PSS:
      {
	RT F = c1c2 * (CGAL_NTS sqrt(D1D2) + b1b2 - a1a2);
	return I + p_ref().y() * hw() - RT(int(_r)) * CGAL_NTS sqrt(F);
      }
      break;
    case SSS:
      {
	return ( cy[0] * CGAL_NTS sqrt(D[0]) +
		 cy[1] * CGAL_NTS sqrt(D[1]) +
		 cy[2] * CGAL_NTS sqrt(D[2]) );
      }
      break;
    default:
      return RT(0);
    }
  }

  RT hw() const {
    switch (v_type) {
    case PPP:
      {
	return uz;
      }
      break;
    case PPS:
      {
	return X;
      }
      break;
    case PSS:
      {
	return CGAL_NTS sqrt(D1D2) - a1a2 - b1b2;
      }
      break;
    case SSS:
      {
	return ( cz[0] * CGAL_NTS sqrt(D[0]) +
		 cz[1] * CGAL_NTS sqrt(D[1]) +
		 cz[2] * CGAL_NTS sqrt(D[2]) );
      }
      break;
    default:
      return RT(1);
    }
  }

  FT squared_radius() const {
    switch (v_type) {
    case PPP:    case PPS:    case PSS:
      {
	Point p_ref;
	if ( p_.is_point() ) {
	  p_ref = p_.point();
	} else if ( q_.is_point() ) {
	  p_ref = q_.point();
	} else {
	  CGAL_assertion( r_.is_point() );
	  p_ref = r_.point();
	}

	FT dx2 = CGAL_NTS square(x() - p_ref.x());
	FT dy2 = CGAL_NTS square(y() - p_ref.y());
	return dx2 + dy2;
      }
      break;
    case SSS:
      {
#if 1
	RT a, b, c;
	compute_supporting_line(p_.segment(), a, b, c);
	Line l(a, b, c);
	Homogeneous_point_2 q = compute_projection(l, point());
#else
	Line l = p_.segment().supporting_line();
	Point q = l.projection(point());
#endif
	FT dx2 = CGAL_NTS square(x() - q.x());
	FT dy2 = CGAL_NTS square(y() - q.y());
	return dx2 + dy2;
      }
      break;
    default:
      return FT(0);
    }
  }


  inline Point  point() const {
    if ( is_degenerate_Voronoi_circle() ) {
      return p_ref();
    }
    
    return Point(x(), y());
  }

  inline Circle circle() const {
    return Circle(point(), squared_radius());
  }

  inline vertex_t type() const { return v_type; }

public:
  Svd_voronoi_vertex_2(const Site& p, const Site& q, const Site &r)
    : p_(p), q_(q), r_(r)
  {
    compute_vertex(p, q, r);
  }

  //--------------------------------------------------------------------------

  template<class Method_tag>
  Sign incircle(const Site& t, Method_tag tag) const 
  {
    if ( t.is_point() ) {
      return incircle(t.point(), tag);
    }
    return incircle(t.segment(), tag);
  }

  //--------------------------------------------------------------------------


  template<class Method_tag>
  Orientation orientation(const Line& l, Method_tag tag) const 
  {
    Orientation o(COLLINEAR);
    switch ( v_type ) {
    case PPP:
      o = orientation(l, PPP_Type(), tag);
      break;
    case PPS:
      o = orientation(l, PPS_Type(), tag);
      break;
    case PSS:
      o = orientation(l, PSS_Type(), tag);
      break;
    case SSS:
      o = orientation(l, SSS_Type(), tag);
      break;
    }

    return o;
  }

  template<class Method_tag>
  Oriented_side oriented_side(const Line& l, Method_tag tag) const
  {
#if 0
    return l.oriented_side(point());
#else
    Orientation o = orientation(l, tag);

    if ( o == COLLINEAR ) { return ON_ORIENTED_BOUNDARY; }
    return ( o == LEFT_TURN ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
#endif
  }

  //--------------------------------------------------------------------------

private:
  const Site& p_, q_, r_;

  vertex_t v_type;

  // index that indicates the refence point for the case PPS
  short pps_idx;

  // the case ppp
  RT ux, uy, uz;

  // the case pss
  RT I, J, c1c2, a1a2, b1b2, D1D2;
  Sign _s, _r;

  // the case pps
  RT x_, y_, X, S; // also I, J from case "ssp"

  // the case sss
  RT cx[3], cy[3], cz[3], D[3];


};





//***********************************************************************
//***********************************************************************
//                            CONSTRUCTIONS
//***********************************************************************
//***********************************************************************


//-----------------------------------------------------------------------
//                  Segment Voronoi diagram vertex
//-----------------------------------------------------------------------

template < class K >
class Construct_svd_vertex_2
{
public:
  typedef typename K::Site_2        Site;
  typedef Svd_voronoi_vertex_2<K>   Voronoi_vertex;
  typedef typename K::Point_2       Point;

  inline Point operator() ( const Site& s1,
			    const Site& s2,
			    const Site& s3) const
  {
    Voronoi_vertex v(s1, s2, s3);
    return v.point();
  }
};


//-----------------------------------------------------------------------
//                  Segment Voronoi diagram circle
//-----------------------------------------------------------------------


template < class Gt >
class Construct_svd_circle_2
{
public:
  typedef typename Gt::Site_2                 Site;
  typedef Svd_voronoi_vertex_2<Gt>            Voronoi_vertex;
  typedef typename Gt::Circle_2               Circle;

  inline Circle operator() ( const Site& s1,
			     const Site& s2,
			     const Site& s3) const
  {
    Voronoi_vertex v(s1, s2, s3);
    return v.circle();
  }
};



//-----------------------------------------------------------------------
//                    Segment Voronoi diagram bisector
//-----------------------------------------------------------------------


template < class Gt >
class Construct_svd_bisector_2
{
public:
  typedef typename Gt::Site_2        Site;
  typedef typename Gt::Point_2       Point;
  typedef typename Gt::Line_2        Line;

  inline Line operator() ( const Site& p, const Site& q) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    if ( p.is_point() && q.is_point() ) {
      Point mid = midpoint(p.point(), q.point());
      Line l(p.point(), q.point());
      return l.perpendicular(mid);
    }
    if ( p.is_segment() && q.is_point() ) {
      // in this case q has to be one of the two endpoints of the
      // segment p...
      Line l = p.segment().supporting_line();
      return l.perpendicular(q.point());
    }
    // in this case p has to be one of the two endpoints of the
    // segment q...
    Line l = q.segment().supporting_line();
    return l.perpendicular(p.point());
  }
};

//-----------------------------------------------------------------------
//                 Segment Voronoi diagram bisector ray
//-----------------------------------------------------------------------

template < class Gt >
class Construct_svd_bisector_ray_2
{
public:
  typedef typename Gt::Site_2                 Site;
  typedef typename Gt::Point_2                Point;
  typedef typename Gt::Line_2                 Line;
  typedef typename Gt::Ray_2                  Ray;
  typedef CGAL::Construct_svd_vertex_2<Gt>    Construct_svd_vertex_2;

  inline Ray operator() ( const Site& p, const Site& q,
			  const Site& r) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    Point v = Construct_svd_vertex_2()(p, q, r);
    Point p1, p2;
    if ( p.is_point() && q.is_point() ) {
      p1 = q.point();
      p2 = p.point();
    } else if ( p.is_point() && q.is_segment() ) {
      CGAL_assertion( p.point() == q.source() ||
		      p.point() == q.target() );
      p1 = (p.point() == q.source()) ? q.target() : q.source();
      p2 = p.point();
    } else {
      // p is a segment and q a point
      p1 = q.point();
      p2 = (q.point() == p.source()) ? p.target() : p.source();
    }
    Line l(p1, p2);
    //    Point base = p.is_point() ? p.point() : q.point();
    Line lperp = l.perpendicular( v );
    return Ray(v, lperp.direction());
  }
};


//-----------------------------------------------------------------------
//              Segment Voronoi diagram bisector segment
//-----------------------------------------------------------------------


template < class Gt >
class Construct_svd_bisector_segment_2
{
public:
  typedef typename Gt::Site_2                     Site;
  typedef typename Gt::Point_2                    Point;
  typedef typename Gt::Line_2                     Line;
  typedef typename Gt::Ray_2                      Ray;
  typedef typename Gt::Segment_2                  Segment;
  typedef CGAL::Parabola_segment_2<Gt>            Parabola_segment;

  typedef CGAL::Construct_svd_vertex_2<Gt>        Construct_svd_vertex_2;
  typedef typename Gt::Are_same_points_2          Are_same_points_2;

  inline Object operator() ( const Site& p, const Site& q,
			     const Site& r, const Site& s) const
  {
    Construct_svd_vertex_2 circumcenter;
    Point vpqr = circumcenter(p, q, r);
    Point vqps = circumcenter(q, p, s);


    Are_same_points_2 are_same_points;

    if ( (p.is_point() && q.is_point()) ||
	 (p.is_segment() && q.is_segment()) ) {
      Segment vorseg(vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    if ( p.is_point() ) {
      // check is p is an endpoint of q
      if (  are_same_points( p.point(), q.segment().source() ) ||
	    are_same_points( p.point(), q.segment().target() )  ) {
	Segment vorseg(vpqr, vqps);
	return CGAL::make_object(vorseg);
      }
      Line l = q.segment().supporting_line();
      Parabola_segment vorseg(p.point(), l, vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    // check is q is an endpoint of p
    if ( q.point() == p.segment().source() ||
	 q.point() == p.segment().target() ) {
      Segment vorseg(vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    Line l = p.segment().supporting_line();
    Parabola_segment vorseg(q.point(), l, vpqr, vqps);
    return CGAL::make_object(vorseg);
  }
};

//-----------------------------------------------------------------------


CGAL_END_NAMESPACE



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_CONSTRUCTIONS_C2_H
