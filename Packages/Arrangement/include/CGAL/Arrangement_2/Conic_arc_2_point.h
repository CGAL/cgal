// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Arrangement_2/Conic_arc_2_point.h
// package       : Arrangement (2.62)
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_CONIC_ARC_2_POINT_H
#define CGAL_CONIC_ARC_2_POINT_H

CGAL_BEGIN_NAMESPACE

template <class R>
class Point_2_ex : public Point_2<R>
{
 public:

   typedef typename R::FT      NT;
 
  enum Type
  {
    User_defined,
    Tangency,
    Intersection_exact,
    Intersection_approx,
    Ray_shooting_exact,
    Ray_shooting_approx
  };

 private:

  Type            _type;
  int             conic_id1;
  int             conic_id2;
  APNT            x_app;
  APNT            x_err;
  int             x_deg;
  std::vector<NT> x_coeffs;
  APNT            y_app;
  APNT            y_err;
  int             y_deg;
  std::vector<NT> y_coeffs;

 public:

  // Constructors.
  Point_2_ex () :
    CGAL::Point_2<R>(),
    _type(User_defined),
    conic_id1(0),
    conic_id2(0),
    x_app(0),
    x_err(0),
    x_deg(0),
    y_app(0),
    y_err(0),
    y_deg(0)
  {}

  Point_2_ex (const NT& hx, const NT& hy, const NT& hz,
	      const Type& type) :
    CGAL::Point_2<R>(hx,hy,hz),
    _type(type),
    conic_id1(0),
    conic_id2(0),
    x_deg(0),
    y_deg(0)
  {
    x_app = TO_APNT(x());
    x_err = 2 * x().get_double_error();
    y_app = TO_APNT(y());
    y_err = 2 * y().get_double_error();
  }

  Point_2_ex (const NT& hx, const NT& hy,
	      const Type& type,
	      const int& id1 = 0, const int& id2 = 0) :
    CGAL::Point_2<R>(hx,hy),
    _type(type),
    conic_id1(id1),
    conic_id2(id2),
    x_deg(0),
    y_deg(0)
  {
    x_app = TO_APNT(x());
    x_err = 2 * x().get_double_error();
    y_app = TO_APNT(y());
    y_err = 2 * y().get_double_error();
  }

  Point_2_ex (const NT& hx, const NT& hy) :
    CGAL::Point_2<R>(hx,hy),
    _type(User_defined),
    conic_id1(0),
    conic_id2(0),
    x_deg(0),
    y_deg(0)
  {
    x_app = TO_APNT(x());
    x_err = 2 * x().get_double_error();
    y_app = TO_APNT(y());
    y_err = 2 * y().get_double_error();
  }

  // Attach the generating polynomials information.
  void attach_polynomials (const int& _x_deg, std::vector<NT> _x_coeffs,
			   const int& _y_deg, std::vector<NT> _y_coeffs)
  {
    // Copy the polynomials.
    x_deg = _x_deg;
    x_coeffs = _x_coeffs;
    y_deg = _y_deg;
    y_coeffs = _y_coeffs;

    // Adjust the degrees.
    if (x_deg > 0)
    {
      while (x_deg > 0 && x_coeffs[x_deg] == 0)
	x_deg--;
      CGAL_assertion(x_deg > 0);
    }

    if (y_deg > 0)
    {
      while (y_deg > 0 && y_coeffs[y_deg] == 0)
	y_deg--;
      CGAL_assertion(y_deg > 0);
    }

    // Adjust the error bounds for roots of polynomials of degree > 2.
    double  slope;
    NT      val;
    int     k;

    if (x_deg > 2)
    {
      // Calculate p'(x).
      slope = x_deg*TO_APNT(x_coeffs[x_deg]);
      for (k = x_deg - 1; k >= 1 ; k--)
	slope = slope*x_app + k*TO_APNT(x_coeffs[k]);

      slope = APNT_ABS(slope);
      if (slope > 1)
	slope = 1;

      if (slope > APNT_EPS)
      {
	// Calculate p(x_app).
	val = x_coeffs[x_deg];
	for (k = x_deg - 1; k >= 0; k--)
	  val = val*x() + x_coeffs[k];

	// Now since |x - x_app|*p'(x) ~= p(x_app) we get:
	x_err = TO_APNT(x_deg*val / NT(slope));
      }
    }

    if (y_deg > 2)
    {
      // Calculate p'(y).
      slope = x_deg*TO_APNT(y_coeffs[y_deg]);
      for (k = y_deg - 1; k >= 1 ; k--)
	slope = slope*y_app + k*TO_APNT(y_coeffs[k]);

      slope = APNT_ABS(slope);
      if (slope > 1)
	slope = 1;

      if (slope > APNT_EPS)
      {
	// Calculate p(y_app).
	val = y_coeffs[y_deg];
	for (k = y_deg - 1; k >= 0; k--)
	  val = val*y() + y_coeffs[k];

	// Now since |y - y_app|*p'(y) ~= p(y_app) we get:
	y_err = TO_APNT(y_deg*val / NT(slope));
      }
    }

    return;
  }

  // Get the approximated co-ordinates.
  const APNT& x_approx () const
  {
    return (x_app);
  }

  const APNT& y_approx () const
  {
    return (y_app);
  }

  // Get the approximation error of the co-ordinates.
  const APNT& x_approx_err () const
  {
    return (x_err);
  }

  const APNT& y_approx_err () const
  {
    return (y_err);
  }

  // Check if the point is approximate.
  bool is_approximate () const
  {
    return (_type == Intersection_approx ||
	    _type == Ray_shooting_approx);
  }

  // Check if the given conic generates the points.
  bool is_generating_conic_id (const int& id) const
  {
    return (id == conic_id1 || id == conic_id2);
  }

  // Compare the co-ordinates of two given points.
  Comparison_result compare_x (const Point_2_ex<R>& p) const
  {
#ifdef CGAL_CONIC_ARC_USE_FILTER
    if (APNT_ABS(x_app - p.x_app) >= (x_err + p.x_err))
    {
      // Using the approximated values is safe in this case!
      if (x_app > p.x_app)
	return (LARGER);
      else if (x_app < p.x_app)
	return (SMALLER);
    }
#endif

    if ((is_approximate() || p.is_approximate()) &&
	eps_compare<APNT>(TO_APNT(x()), TO_APNT(p.x())) == EQUAL)
    {
      if (x_deg > 0 && (p.x_deg == 0 || !p.is_approximate()))
      {
	if (_is_root (p.x(), x_coeffs, x_deg))
	  return (EQUAL);
      }
      else if ((x_deg == 0 || ! is_approximate()) && p.x_deg > 0)
      {
	if (_is_root (x(), p.x_coeffs, p.x_deg))
	  return (EQUAL);
      }
      else if (x_deg > 0 && p.x_deg > 0)
      {
	if ((conic_id1 == p.conic_id1 && conic_id2 == p.conic_id2) ||
	    (conic_id1 == p.conic_id2 && conic_id2 == p.conic_id1))
	  return (EQUAL);

	std::vector<NT> gcd_coeffs;
	int             deg_gcd;

	gcd_polynomials<NT> (x_coeffs, x_deg,
			     p.x_coeffs, p.x_deg,
			     gcd_coeffs, deg_gcd);

	if (deg_gcd == x_deg)
	  return (EQUAL);

	if (deg_gcd > 0 && _is_approx_root (x_app, x_err, 
					    gcd_coeffs, deg_gcd))
	{
	  return (EQUAL);
	}

	// The two roots are very close: perform exhaustive comparison.
	return (compare_roots<NT>(x_coeffs, x_deg, 
				  NT(x_app), NT(x_err),
				  p.x_coeffs, p.x_deg, 
				  NT(p.x_app), NT(p.x_err)));
      }
    }

    return (CGAL::compare (x(), p.x()));
  }

  Comparison_result compare_y (const Point_2_ex<R>& p) const
  {
#ifdef CGAL_CONIC_ARC_USE_FILTER
    if (APNT_ABS(y_app - p.y_app) >= (y_err + p.y_err))
    {
      // Using the approximated values is safe in this case!
      if (y_app > p.y_app)
	return (LARGER);
      else if (y_app < p.y_app)
	return (SMALLER);
    }
#endif

    if ((is_approximate() || p.is_approximate()) &&
	eps_compare<APNT>(TO_APNT(y()), TO_APNT(p.y())) == EQUAL)
    {
      if (y_deg > 0 && (p.y_deg == 0 || ! p.is_approximate()))
      {
	if (_is_root (p.y(), y_coeffs, y_deg))
	  return (EQUAL);
      }
      else if ((y_deg == 0 || ! is_approximate()) && p.y_deg > 0)
      {
	if (_is_root (y(), p.y_coeffs, p.y_deg))
	  return (EQUAL);
      }
      else if (y_deg > 0 && p.y_deg > 0)
      {
	if ((conic_id1 == p.conic_id1 && conic_id2 == p.conic_id2) ||
	    (conic_id1 == p.conic_id2 && conic_id2 == p.conic_id1))
	  return (EQUAL);

	std::vector<NT> gcd_coeffs;
	int             deg_gcd;

	gcd_polynomials<NT> (y_coeffs, y_deg,
			     p.y_coeffs, p.y_deg,
			     gcd_coeffs, deg_gcd);

	if (deg_gcd == y_deg)
	  return (EQUAL);

	if (deg_gcd > 0 && _is_approx_root (y_app, y_err, 
					    gcd_coeffs, deg_gcd))
	{
	  return (EQUAL);
	}
	
	// The two roots are very close: perform exhaustive comparison.
	return (compare_roots<NT>(y_coeffs, y_deg, 
				  NT(y_app), NT(y_err),
				  p.y_coeffs, p.y_deg, 
				  NT(p.y_app), NT(p.y_err)));
      }
    }

    return (CGAL::compare (y(), p.y()));
  }
  
  // Equality operators.
  bool equals (const Point_2_ex<R>& p) const
  {
    return (compare_x(p) == EQUAL && compare_y(p) == EQUAL);
  }

  bool operator== (const Point_2_ex<R>& p) const
  {
    return (compare_x(p) == EQUAL && compare_y(p) == EQUAL);
  }

  bool operator!= (const Point_2_ex<R>& p) const
  {
    return (compare_x(p) != EQUAL || compare_y(p) != EQUAL);
  }

  // Copmare two points lexicographically.
  Comparison_result compare_lex_xy (const Point_2_ex<R>& p) const
  {
    Comparison_result   res = this->compare_x (p);

    if (res != EQUAL)
      return (res);
    
    return (this->compare_y (p));
  }

  Comparison_result compare_lex_yx (const Point_2_ex<R>& p) const
  {
    Comparison_result   res = this->compare_y (p);

    if (res != EQUAL)
      return (res);
    
    return (this->compare_x (p));
  }

  // Reflect a point.
  Point_2_ex<R> reflect_in_y () const
  {
    Point_2_ex<R> ref_point (-hx(), hy(), hw(), _type);
    int           i;

    if (conic_id1 != 0)
      ref_point.conic_id1 = conic_id1 ^ REFLECT_IN_Y;
    if (conic_id2 != 0)
      ref_point.conic_id2 = conic_id2 ^ REFLECT_IN_Y;

    ref_point.x_app = -x_app;
    ref_point.x_err = x_err;
    ref_point.x_deg = x_deg;
    if (x_deg > 0)
    {
      ref_point.x_coeffs.resize (x_deg + 1);
      for (i = 0; i <= x_deg; i++)
	ref_point.x_coeffs[i] = 
	  (i % 2 == 0) ? x_coeffs[i] : -x_coeffs[i];
    }

    ref_point.y_app = y_app;
    ref_point.y_err = y_err;
    ref_point.y_deg = y_deg;
    ref_point.y_coeffs = y_coeffs;

    return (ref_point);
  }

  Point_2_ex<R> reflect_in_x_and_y () const
  {
    Point_2_ex<R> ref_point (-hx(), -hy(), hw(), _type);
    int           i;

    if (conic_id1 != 0)
      ref_point.conic_id1 = conic_id1 ^ (REFLECT_IN_X | REFLECT_IN_Y);
    if (conic_id2 != 0)
      ref_point.conic_id2 = conic_id2 ^ (REFLECT_IN_X | REFLECT_IN_Y);

    ref_point.x_app = -x_app;
    ref_point.x_err = x_err;
    ref_point.x_deg = x_deg;
    if (x_deg > 0)
    {
      ref_point.x_coeffs.resize(x_deg + 1);
      for (i = 0; i <= x_deg; i++)
	ref_point.x_coeffs[i] = 
	  (i % 2 == 0) ? x_coeffs[i] : -x_coeffs[i];
    }

    ref_point.y_app = -y_app;
    ref_point.y_err = y_err;
    ref_point.y_deg = y_deg;
    if (y_deg > 0)
    {
      ref_point.y_coeffs.resize(y_deg + 1);
      for (i = 0; i <= y_deg; i++)
	ref_point.y_coeffs[i] = 
	  (i % 2 == 0) ? y_coeffs[i] : -y_coeffs[i];
    }

    return (ref_point);
  }

private:

  bool _is_root (const NT& z, 
		 const std::vector<NT>& p, const int& deg) const
  {
    NT   val = p[deg];
    int  k;

    for (k = deg - 1; k >= 0 ; k--)
      val = val*z + p[k];

    return (CGAL::compare (val, NT(0)) == EQUAL);
  }

  bool _is_approx_root (const APNT& z_app, const APNT& z_err, 
			const std::vector<NT>& p, const int& deg) const
  {
    APNT   val;
    APNT   slope;
    int    k;

    // Compute p(z).
    val = TO_APNT(p[deg]);
    for (k = deg - 1; k >= 0 ; k--)
      val = val*z_app + TO_APNT(p[k]);

    // Compute p'(z).
    slope = deg*TO_APNT(p[deg]);
    for (k = deg - 1; k >= 1 ; k--)
      slope = slope*z_app + k*TO_APNT(p[k]);
    slope = APNT_ABS(slope);

    // z approximates a root of p if p(z) < err(z)*p'(z).
    return (APNT_ABS(val) < z_err*slope);
  }

};

CGAL_END_NAMESPACE

#endif

