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
// file          : include/CGAL/Arrangement_2/Conic_arc_2_eq.h
// package       : Arrangement (2.62)
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_CONIC_ARC_2_EQ_H
#define CGAL_CONIC_ARC_2_EQ_H

#include <CGAL/double.h>
#include <CGAL/Arrangement_2/Sturm_seq.h>

// if we use a LEDA version without namespaces
// we have to define a few macros
#if !defined(LEDA_NAMESPACE)
#define LEDA_BEGIN_NAMESPACE
#define LEDA_END_NAMESPACE
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE

typedef double      APNT;
#define APNT_EPS    0.0000001
#define COMP_EPS    0.0001
#define TO_APNT(x)  CGAL::to_double(x)
#define APNT_ABS(x) fabs(x)

// ----------------------------------------------------------------------------
// Compare two approximate values.
//
template <class APNT>
Comparison_result eps_compare (const APNT& val1, const APNT& val2,
			       const APNT& _eps = APNT_EPS)
{
  APNT denom = APNT_ABS(val1 - val2);
  APNT numer = APNT_ABS(val1) + APNT_ABS(val2);

  if (numer < CGAL::sqrt(_eps))
  {
    // The two numbers are very close to zero:
    if (denom < _eps)
      return (EQUAL);
    else
      return ((val1 < val2) ? SMALLER : LARGER);
  }
  else
  {
    //                  | x - y |
    //  x = y    <==>  ----------- < _eps
    //                  |x| + |y|
    if (denom / numer < _eps)
      return (EQUAL);
    else
      return ((val1 < val2) ? SMALLER : LARGER);
  }
}

/*****************************************************************************/
// Some functions for quick-solving equations when working with doubles.

//--------------------------------------------------------------------
// Solve the cubic equation: a*x^3 + b*x^2 + c*x + d = 0 (a != 0),
// with double precision floats and using Tartaglia's method.
// The equation is assumed not to have multiple roots.
// The function returns the number of real solutions.
// The roots area must be at least of size 3.
//
template <class NT>
int _Tartaglia_cubic_eq (const NT& a, 
			 const NT& b, 
			 const NT& c, 
			 const NT& d,
			 NT *roots)
{
  // Some constants.
  static const NT    _Zero = 0;
  static const NT    _One = 1;
  static const NT    _Two = 2;
  static const NT    _Three = 3;
  static const NT    _Nine = 9;
  static const NT    _Twenty_seven = 27;
  static const NT    _Fifty_four = 54;
  static const NT    _One_third = _One / _Three;

  // Normalize the equation to obtain: x^3 + a2*x^2 + a1*x^1 + a0 = 0.
  const NT    a2 = b / a;
  const NT    a1 = c / a;
  const NT    a0 = d / a;

  // Caclculate the "cubic discriminant" D.
  const NT    Q = (_Three*a1 - a2*a2) / 
                      _Nine;
  const NT    R = (_Nine*a1*a2 - _Twenty_seven*a0 - _Two*a2*a2*a2) / 
                      _Fifty_four;
  const NT    D = Q*Q*Q + R*R;
  NT          S, T;

  if (D > _Zero)
  {
    // One real-valued root (the other two are complex conjugates):
    NT       temp = R + CGAL::sqrt(D);

    if (temp > _Zero)
      S = ::pow(temp, _One_third);
    else
      S = - ::pow(-temp, _One_third);

    temp = R - CGAL::sqrt(D);
    if (temp > _Zero)
      T = ::pow(temp, _One_third);
    else
      T = - ::pow(-temp, _One_third);
    
    roots[0] = S + T - a2/_Three;
    return (1);
  }
  else
  {
    // We have three real-valued roots:
    const NT  phi = ::atan2 (CGAL::sqrt(-D), R) / _Three;
    const NT  sin_phi = ::sin(phi);
    const NT  cos_phi = ::cos(phi);

    roots[0] = _Two*CGAL::sqrt(-Q)*cos_phi - a2/_Three;
    roots[1] = -CGAL::sqrt(-Q)*cos_phi - 
               CGAL::sqrt(-_Three*Q)*sin_phi - a2/_Three;
    roots[2] = -CGAL::sqrt(-Q)*cos_phi + 
               CGAL::sqrt(-_Three*Q)*sin_phi - a2/_Three;

    return (3);
  }
}

//--------------------------------------------------------------------
// Calculate the square root of a complex number.
//
template <class NT>
void _complex_sqrt (const NT& real, const NT& imag,
		    NT& sqrt_re, NT& sqrt_im,
		    const NT& _eps)
{
  static const NT    _Zero = 0;
  static const NT    _Two = 2;

  // Check the special case where imag=0:
  if (eps_compare<NT>(imag, _Zero, _eps) == EQUAL)
  {
    if (real >= _Zero)
    {
      sqrt_re = CGAL::sqrt(real);
      sqrt_im = _Zero;
    }
    else
    {
      sqrt_re = _Zero;
      sqrt_im = CGAL::sqrt(-real);
    }

    return;
  }

  // We seek x,y such that: (x^2 - y^2 = 0) and (2xy = imag).
  // So, we get the equation: x^4 - real*x^2 - (imag/2)^2 = 0.
  // This equation has one positive root (and one irrelevant negative root).
  const NT  sqrt_disc = CGAL::sqrt(real*real + imag*imag);
  const NT  x_2 = (real + sqrt_disc) / _Two;
  
  // Calculate on of the roots (the other root is -sqrt_re - i*sqrt_im).
  sqrt_re = CGAL::sqrt(x_2);
  sqrt_im = imag / (_Two*sqrt_re);

  return;
}

//--------------------------------------------------------------------
// Solve the quartic equation: a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 (a != 0),
// with double precision floats and using Ferrari's method.
// The equation is assumed not to have multiple roots.
// The function returns the number of real solutions.
// The roots area must be at least of size 4.
//
template <class NT>
int _Ferrari_quartic_eq (const NT& a, 
			 const NT& b, 
			 const NT& c,
			 const NT& d, 
			 const NT& e,
			 NT *roots,
			 const NT& _eps)
{
  // Some constants.
  static const NT    _Zero = 0;
  static const NT    _One = 1;
  static const NT    _Two = 2;
  static const NT    _Three = 3;
  static const NT    _Four = 4;
  static const NT    _Eight = 8;
  static const NT    _One_quarter = _One/_Four;
  static const NT    _One_half = _One/_Two;
  static const NT    _Three_quarters = _Three/_Four;

  // Normalize the equation to obtain: x^4 + a3*x^3 + a2*x^2 + a1*x^1 + a0 = 0.
  const NT    a3 = b / a;
  const NT    a2 = c / a;
  const NT    a1 = d / a;
  const NT    a0 = e / a;

  // Find a real root of the resolvent cubic equation:
  // y^3 - a2*y^2 + (a1*a3 - 4*a0)*y + (2*a2*a0 - a1^2 - a0*a3^2) = 0.
  NT          cubic_roots[3];
  NT          y1;

  _Tartaglia_cubic_eq (_One, 
		       -a2, 
		       a1*a3 - _Four*a0, 
		       _Four*a2*a0 - a1*a1 - a3*a3*a0,
		       cubic_roots);
  y1 = cubic_roots[0];

  // Calculate R^2 and act accordingly.
  const NT    R_square = _One_quarter*a3*a3 - a2 + y1;
  int             n_roots = 0;

  if (eps_compare<NT>(R_square, _Zero, _eps) == EQUAL)
  {
    // R is _Zero: In this case, D and E must be real-valued.
    const NT    disc_square = y1*y1 - _Four*a0;
      
    if (disc_square < _Zero)
    {
      // No real roots at all.
      return (0);
    }

    const NT    disc = CGAL::sqrt(disc_square);
    const NT    D_square = _Three_quarters*a3*a3 + _Two*(disc - a2);
    const NT    E_square = _Three_quarters*a3*a3 - _Two*(disc + a2);

    if (D_square > _Zero)
    {
      // Add two real-valued roots.
      const NT    D = CGAL::sqrt(D_square);
	  
      roots[n_roots] = -_One_quarter*a3 + _One_half*D;
      n_roots++;
      roots[n_roots] = -_One_quarter*a3 - _One_half*D;
      n_roots++;
    }

    if (E_square > _Zero)
    {
      // Add two real-valued roots.
      const NT    E = CGAL::sqrt(E_square);
      
      roots[n_roots] = -_One_quarter*a3 + _One_half*E;
      n_roots++;
      roots[n_roots] = -_One_quarter*a3 - _One_half*E;
      n_roots++;
    }
  }
  else if (R_square > _Zero)
  {
    // R is real-valued: In this case, D and E must also be real-valued.
    const NT    R = CGAL::sqrt(R_square);
    const NT    D_square = _Three_quarters*a3*a3 - R_square - _Two*a2 +
      (_Four*a3*a2 - _Eight*a1 - a3*a3*a3) / (_Four*R);
    const NT    E_square = _Three_quarters*a3*a3 - R_square - _Two*a2 -
      (_Four*a3*a2 - _Eight*a1 - a3*a3*a3) / (_Four*R);

    if (D_square > _Zero)
    {
      // Add two real-valued roots.
      const NT    D = CGAL::sqrt(D_square);
      
      roots[n_roots] = -_One_quarter*a3 + _One_half*(R + D);
      n_roots++;
      roots[n_roots] = -_One_quarter*a3 + _One_half*(R - D);
      n_roots++;
    }

    if (E_square > _Zero)
    {
      // Add two real-valued roots.
      const NT    E = CGAL::sqrt(E_square);
      
      roots[n_roots] = -_One_quarter*a3 + _One_half*(E - R);
      n_roots++;
      roots[n_roots] = -_One_quarter*a3 - _One_half*(E + R);
      n_roots++;
    }
  }
  else
  {
      // R is imaginary: Calclualte D and E as complex numbers.
      const NT    R_im = CGAL::sqrt(-R_square);
      const NT    R_inv_im = -_One / R_im;      // The inverse of R.
      const NT    D_square_re = _Three_quarters*a3*a3 - R_square - _Two*a2;
      const NT    D_square_im = 
	(a3*a2 - _Two*a1 - _One_quarter*a3*a3*a3) * R_inv_im;
      NT          D_re, D_im;
      const NT    E_square_re = D_square_re;
      const NT    E_square_im = -D_square_im;
      NT          E_re, E_im;

      _complex_sqrt (D_square_re, D_square_im,
		     D_re, D_im, _eps);
      _complex_sqrt (E_square_re, E_square_im,
		     E_re, E_im, _eps);

      // Use only the combinations that yield real-valued numbers.
      if (eps_compare<NT>(R_im, -D_im, _eps) == EQUAL) // (R_im == -D_im)
      {
	roots[n_roots] = -_One_quarter*a3 + _One_half*D_re;
	n_roots++;
      }
      if (eps_compare<NT>(R_im, D_im, _eps) == EQUAL) // (R_im == D_im)
      {
	roots[n_roots] = -_One_quarter*a3 - _One_half*D_re;
	n_roots++;
      }

      if (eps_compare<NT>(R_im, E_im, _eps) == EQUAL) // (R_im == E_im)
      {
	roots[n_roots] = -_One_quarter*a3 + _One_half*E_re;
	n_roots++;
      }
      if (eps_compare<NT>(R_im, -E_im, _eps) == EQUAL) // (R_im == -E_im)
      {
	roots[n_roots] = -_One_quarter*a3 - _One_half*D_re;
	n_roots++;
      }
  }

  return (n_roots);
}

/*****************************************************************************/

// ----------------------------------------------------------------------------
// Solve the quadratic equation: a*x^2 + b*x + c = 0.
// The function returns the number of distinct solutions.
// The roots area must be at least of size 2.
// 
template <class NT>
int solve_quadratic_eq (const NT& a, const NT& b, const NT& c,
			NT* roots, int* mults)
{
  static const NT _zero = 0;
  static const NT _two = 2;
  static const NT _four = 4;

  if (a == _zero)
  {
    // Solve a linear equation.
    if (b != _zero)
    {
      roots[0] = -c/b;
      mults[0] = 1;
      return (1);
    }
    else
      return (0);
  }

  // Act according to the discriminant.
  const NT        disc = b*b - _four*a*c;

  if (disc < _zero)
  {
    // No real roots.
    return (0);
  }
  else if (disc == _zero)
  {
    // One real root with mutliplicity 2.
    roots[0] = roots[1] = -b/(_two*a);
    mults[0] = 2;
    return (1);
  }
  else
  {
    // Two real roots.
    const NT      sqrt_disc = CGAL::sqrt(disc);

    roots[0] = (sqrt_disc - b)/(_two*a);
    mults[0] = 1;
    roots[1] = -(sqrt_disc + b)/(_two*a);
    mults[1] = 1;
    return (2);
  }
}

// ----------------------------------------------------------------------------
// Solve the cubic equation: _a*x^3 + _b*x^2 + _c*x + _d = 0.
// Notice this is an auxiliary function (used only by solve_quartic_eq()).
// The function returns the number of distinct solutions.
// The roots area must be at least of size 3.
// 
template <class NT>
static int _solve_cubic_eq (const NT& _a, const NT& _b,
			    const NT& _c, const NT& _d,
			    NT* roots, int* mults,
			    int& n_approx)
{
  n_approx = 0;

  // In case we have an equation of a lower degree:
  static const NT _zero = 0;

  if (_a == _zero)
  {
    // We have to solve _b*x^2 + _c*x + _d = 0.
    return (solve_quadratic_eq<NT> (_b, _c, _d,
				    roots, mults));
  }

  // Normalize the equation, so that the leading coefficient is 1 and we have:
  //  x^3 + b*x^2 + c*x + d = 0
  //
  const NT  b = _b/_a;
  const NT  c = _c/_a;
  const NT  d = _d/_a;

  // Check whether the equation has multiple roots.
  // If we write: 
  //  p(x) = x^3 + b*x^2 + c*x + d = 0
  //
  // Then:
  //  p'(x) = 3*x^2 + 2*b*x + c
  //
  // We know that there are multiple roots iff GCD(p(x), p'(x)) != 0.
  // In order to check the GCD, let us denote:
  //  p(x) mod p'(x) = A*x + B
  //
  // Then:
  static const NT _two = 2;
  static const NT _three = 3;
  static const NT _nine = 9;
  
  const NT  A = _two*(c - b*b/_three)/_three;
  const NT  B = d - b*c/_nine;

  if (A == _zero)
  {
    // In case A,B == 0, then p'(x) divides p(x).
    // This means the equation has one solution with multiplicity of 3.
    // We can obtain this root using the fact that -b is the sum of p(x)'s 
    // roots.
    if (B == _zero)
    {
      roots[0] = -b / _three;
      mults[0] = 3;
      return (1);
    }
  }
  else
  {
    // Check whether A*x + B divides p'(x).
    const NT   x0 = -B/A;

    if (c + x0*(_two*b + _three*x0) == _zero)
    {
      // Since GCD(p(x), p'(x)) = (x - x0), then x0 is a root of p(x) with
      // multiplicity 2. The other root is obtained from b.
      roots[0] = x0;
      mults[0] = 2;
      roots[1] = -(b + 2*x0);
      mults[1] = 1;
      return (2);
    }
  }
      
  // If we reached here, we can be sure that there are no multiple roots.
  // We use Ferrari's method to approximate the solutions.
  APNT    app_roots[4];
 
  n_approx = _Tartaglia_cubic_eq<APNT> (1,
					TO_APNT(b),
					TO_APNT(c),
					TO_APNT(d),
					app_roots);

  // Copy the approximations.
  for (int i = 0; i < n_approx; i++)
  {
    roots[i] = NT(app_roots[i]);
    mults[i] = 1;
  }

  return (n_approx);
}

// ----------------------------------------------------------------------------
// Solve a quartic equation: _a*x^4 + _b*x^3 + _c*x^2 + _d*x + _e = 0.
// The function returns the number of distinct solutions.
// The roots area must be at least of size 4.
// 
template <class NT>
int solve_quartic_eq (const NT& _a, const NT& _b, const NT& _c, 
		      const NT& _d, const NT& _e,
		      NT* roots, int* mults,
		      int& n_approx)
{
  // Right now we have no approximate solutions.
  n_approx = 0;

  // First check whether we have roots that are 0.
  static const NT _zero = 0;

  if (_e == _zero)
  {
    if (_d == _zero)
    {
      if (_c == _zero)
      {
	if (_b == _zero)
	{
	  // Unless the equation is trivial, we have the root 0,
	  // with multiplicity of 4.
	  if (_a == _zero)
	    return (0);
	  
	  roots[0] = _zero;
	  mults[0] = 4;
	  return (1);
	}
	else
	{
	  // Add the solution 0 (with multiplicity 3), and add the solution
	  // to _a*x + _b = 0.
	  if (_a == _zero)
	    return (0);

	  roots[0] = _zero;
	  mults[0] = 3;
	  roots[1] = -(_b/_a);
	  mults[1] = 1;
	  return (2);
	}
      }
      else
      {
	// Add the solution 0 (with multiplicity 2),
	// and solve _a*x^2 + _b*x + _c*x = 0.
	roots[0] = _zero;
	mults[0] = 2;

	return (1 + solve_quadratic_eq<NT> (_a, _b, _c,
					    roots + 1, mults + 1));
      }
    }
    else
    {
      // Add the solution 0, and solve _a*x^3 + _b*x^2 + _c*x + _d = 0.
      roots[0] = _zero;
      mults[0] = 1;

      return (1 + _solve_cubic_eq<NT> (_a, _b, _c, _d,
				       roots + 1, mults + 1,
				       n_approx));
    }
  }

  // In case we have an equation of a lower degree:
  if (_a == _zero)
  {
    if (_b == _zero)
    {
      // We have to solve _c*x^2 + _d*x + _e = 0.
      return (solve_quadratic_eq<NT> (_c, _d, _e,
				      roots, mults));
    }
    else
    {
      // We have to solve _b*x^3 + _c*x^2 + _d*x + _e = 0.
      return (_solve_cubic_eq<NT> (_b, _c, _d, _e,
				   roots, mults,
				   n_approx));
    }
  }

  // In case we have an equation of the form:
  //  _a*x^4 + _c*x^2 + _e = 0
  //
  // Then by substituting y = x^2, we obtain a quadratic equation:
  //  _a*y^2 + _c*y + _e = 0
  //
  if (_b == _zero && _d == _zero)
  {
    // Solve the equation for y = x^2.
    NT   sq_roots[2];
    int  sq_mults[2];
    int  n_sq;
    int  n = 0;

    n_sq = solve_quadratic_eq<NT> (_a, _c, _e,
				   sq_roots, sq_mults);

    // Convert to roots of the original equation: only positive roots for y
    // are relevant of course.
    for (int j = 0; j < n_sq; j++)
    {
      if (sq_roots[j] > _zero)
      {
	roots[n] = CGAL::sqrt(sq_roots[j]);
	mults[n] = sq_mults[j];
	roots[n+1] = -roots[n];
	mults[n+1] = sq_mults[j];
	n += 2;
      }
    }

    return (n);
  }

  // Normalize the equation, so that the leading coefficient is 1 and we have:
  //  x^4 + b*x^3 + c*x^2 + d*x + e = 0
  // 
  const NT b = _b/_a;
  const NT c = _c/_a;
  const NT d = _d/_a;
  const NT e = _e/_a;

  // Check whether the equation has multiple roots.
  // If we write: 
  //  p(x) = x^4 + b*x^3 + c*x^2 + d*x + e = 0
  //
  // Then:
  //  p'(x) = 4*x^3 + 3*b*x^2 + 2*c*x + d
  //
  // We know that there are multiple roots iff GCD(p(x), p'(x)) != 0.
  // In order to check the GCD, let us denote:
  //  p(x) mod p'(x) = alpha*x^2 + beta*x + gamma
  //
  // Then:
  static const NT _two = 2;
  static const NT _three = 3;
  static const NT _four = 4;
  static const NT _eight = 8;
  static const NT _sixteen = 16;

  const NT alpha = c/_two - _three*b*b/_sixteen;
  const NT beta = _three*d/_four - b*c/_eight;
  const NT gamma = e - b*d/_sixteen;

  if (alpha == _zero)
  {
    if (beta == _zero)
    {
      if (gamma == _zero)
      {
	// p'(x) divides p(x). This can happen only if there is a single
	// root with multiplicity 4.
	// Since b = -(sum of p(x)'s roots) then this root is -b/4.
	roots[0] = -b/_four;
	mults[0] = 4;
	return (1);
      }
    }
    else
    {
      // In case alpha is 0, check whether (beta*x + gamma) divides p'(x) -
      // this happens iff -gamma/beta is a root of p'(x).
      if (beta == _zero)
	return (0);

      const NT x0 = -gamma/beta;

      if (d + x0*(_two*c + x0*(_three*b + x0*_four)) == _zero)
      {
	// We have two ditinct solutions: x0 with multiplicity of 3, and
	// p'(x)/(x - x0)^2 with multiplicity of 1.
	// Since b = -(sum of p(x)'s roots), the other root is -(b + 3*x0).
	roots[0] = x0;
	mults[0] = 3;
	roots[1] = -(b + 3*x0);
	mults[1] = 1;
	return (2);
      }
    }
  }
  else
  {
    // If alpha is not 0, let us denote:
    //  p'(x) mod (alpha*x^2 + beta*x + gamma) = A*x + B
    //
    // And so:
    const NT A = _two*c - (_four*gamma + _three*b*beta)/alpha + 
                 _four*beta*beta/(alpha*alpha);
    const NT B = d - _three*b*gamma/alpha + _four*beta*gamma/(alpha*alpha);

    // In case A,B == 0, then GCD(p(x), p'(x)) = (alpha*x^2 + beta*x + gamma).
    // This means the equation has two distinct solution, each one of them
    // has a multiplicity of 2.
    if (A == _zero)
    {
      if (B == _zero)
      {
	// There are two cases: 
	// - If (alpha*x^2 + beta*x + gamma) has two distinct roots, that 
	//   these roots are roots of p(x), each with multiplicity 2.
	// - If there is one root x0, than this root is a root of p(x) with
	//   multiplicity 3. The other root can be obtained using the fact  
	//   that b = -(sum of p(x)'s roots).	
	int   m  = solve_quadratic_eq<NT> (alpha, beta, gamma,
					   roots, mults);

	if (m == 1)
	{
	  // Add the second root.
	  mults[0] = 3;
	  roots[1] = -(b + 3*roots[0]);
	  mults[1] = 1;
	  m++;
	}
	else if (m == 2)
	{
	  mults[0] = 2;
	  mults[1] = 2;
	}

	return (m);
      }
    }
    else
    {
      // Check whether (A*x + B) divides (alpha*x^2 + beta*x + gamma).
      // This happend iff -B/A is a root of (alpha*x^2 + beta*x + gamma).
      const NT x0 = -B / A;

      if (gamma + x0*(beta + x0*alpha) == _zero)
      {
	// We have one solution with multiplicity of 2.
	roots[0] = x0;
	mults[0] = 2;

	// The other two solutions are obtained from solving:
	// p(x) / (x - x0)^2 = x^2 + (b + 2*x0)*x + (c + 2*b*x0 + 3*x0^2)
	return (1 + 
		solve_quadratic_eq<NT> (NT(1), b + 2*x0, c + x0*(2*b + 3*x0),
					roots + 1, mults + 1));		
      }
    }
  }
  
  // Compute the number of real roots we should compute, using Sturm sequences.
  Polynom<NT>    p;
  p.set (4, _a);
  p.set (3, _b);
  p.set (2, _c);
  p.set (1, _d);
  p.set (0, _e);

  Sturm_seq<NT>  sts (p);
  int            n_expected = sts.sign_changes_at_infinity(-1) -
                              sts.sign_changes_at_infinity(1);
 
  // If we reached here, we can be sure that there are no multiple roots.
  // We use Ferrari's method to approximate the solutions.
  APNT           app_roots[4];
  APNT           _eps = APNT_EPS;

  do
  {
    n_approx = _Ferrari_quartic_eq<APNT> (1,
					  TO_APNT(b),
					  TO_APNT(c),
					  TO_APNT(d),
					  TO_APNT(e),
					  app_roots,
					  _eps);

    _eps *= 10;
  } while (n_approx != n_expected && _eps < 1);
  
  CGAL_assertion(n_approx == n_expected);

  // Copy the approximations.
  for (int i = 0; i < n_approx; i++)
  {
    roots[i] = NT(app_roots[i]);
    mults[i] = 1;
  }
    
  return (n_approx);
}

// ----------------------------------------------------------------------------
// Find the GCD of the two given polynomials.
// 
template <class NT>
void _mod_polynomials (const NT* a_coeffs, const int& a_deg,
		       const NT* b_coeffs, const int& b_deg,
		       NT* r_coeffs, int& r_deg)
{
  const static NT _zero = 0;
  NT              work[5];
  int             deg = a_deg;
  NT              factor;
  int             k;

  // Copy the polynomial a(x) to the work area.
  for (k = 0; k <= a_deg; k++)
    work[k] = a_coeffs[k];

  // Perform long division of a(x) by b(x).
  while (deg >= b_deg)
  {
    factor = work[deg] / b_coeffs[b_deg];

    for (k = b_deg - 1; k >= 0; k--)
      work[deg - b_deg + k] -= factor*b_coeffs[k];

    do
    {
      deg--;
    } while (deg >= 0 && work[deg] == _zero);
  }

  // Copy back the residue.
  for (k = 0; k <= deg; k++)
    r_coeffs[k] = work[k];
  r_deg = deg;

  return;
}

template <class NT>
void gcd_polynomials (const std::vector<NT>& a_coeffs, const int& a_deg,
		      const std::vector<NT>& b_coeffs, const int& b_deg,
		      std::vector<NT>& gcd_coeffs, int& gcd_deg)
{
  // Copy the polynomials to the work areas (make sure that A has a higher
  // degree than B).
  NT  a[5], b[5], r[5];
  NT  *aP = &(a[0]), *bP=&(b[0]), *rP = &(r[0]);
  NT  *swapP;
  int adeg, bdeg, rdeg;
  int k;

  if (a_deg >= b_deg)
  {
    for (k = 0; k <= a_deg; k++)
      aP[k] = a_coeffs[k];
    adeg = a_deg;

    for (k = 0; k <= b_deg; k++)
      bP[k] = b_coeffs[k];
    bdeg = b_deg;
  }
  else
  {
    for (k = 0; k <= b_deg; k++)
      aP[k] = b_coeffs[k];
    adeg = b_deg;

    for (k = 0; k <= a_deg; k++)
      bP[k] = a_coeffs[k];
    bdeg = a_deg;
  }

  // Perform Euclid's algorithm.
  while (true)
  {
    _mod_polynomials (aP, adeg, bP, bdeg,
		      rP, rdeg);

    if (rdeg <= 0)
      break;

    swapP = aP;
    aP = bP;
    adeg = bdeg;
    bP = rP;
    bdeg = rdeg;
    rP = swapP;
  }

  if (rdeg < 0)
  {
    // b(x) divides a(x), so gcd(a,b) = b.
    gcd_coeffs.resize(bdeg + 1);
    for (k = 0; k <= bdeg; k++)
      gcd_coeffs[k] = bP[k];
    gcd_deg = bdeg;
  }
  else
  {
    // The two polynomials have gcd(a,b) = 1.
    gcd_deg = 0;
  }

  return;
}

// ----------------------------------------------------------------------------
// Check whether alpha is a root of the given polynomial p(x).
// If so, compute: q(x) = p(x)/(x - alpha).
// 
template <class NT>
bool factor_root (const NT* p_coeffs, const int& p_deg,
		  const NT& alpha,
		  NT* q_coeffs, int& q_deg)
{
  const static NT _zero = 0;
  NT              coeff = p_coeffs[p_deg];
  int             k;

  for (k = p_deg; k > 0; k--)
  {
    q_coeffs[k-1] = coeff;
    coeff = p_coeffs[k-1] + alpha*coeff;
  }
  
  if (coeff == _zero)
  {
    q_deg = p_deg - 1;
    return (true);
  }
  else
  {
    q_deg = 0;
    return (false);
  }
}

// ----------------------------------------------------------------------------
// Check which one of the polynomial roots is greater than the other.
// 
template <class NT>
Comparison_result compare_roots (const std::vector<NT>& p1_coeffs, 
				 const int& p1_deg,
				 const NT& alpha1, const NT& delta1,  
				 const std::vector<NT>& p2_coeffs,
				 const int& p2_deg,
				 const NT& alpha2, const NT& delta2)
{
  NT           l1 = alpha1 - delta1;
  NT           r1 = alpha1 + delta1;
  NT           l2 = alpha2 - delta2;
  NT           r2 = alpha2 + delta2;

  if (r1 < l2)
    return (SMALLER);
  else if (l1 > r2)
    return (LARGER);

  // Generate Sturm sequences.
  Polynom<NT>    p1;
  Polynom<NT>    p2;
  APNT           max_coeff = 0;
  int            i;

  for (i = p1_deg; i >= 0; i--)
  {
    p1.set (i, p1_coeffs[i]);
    if (TO_APNT(p1_coeffs[i]) > max_coeff)
      max_coeff = TO_APNT(p1_coeffs[i]);
  }
  
  for (i = p2_deg; i >= 0; i--)
  {
    p2.set (i, p2_coeffs[i]);
    if (TO_APNT(p2_coeffs[i]) > max_coeff)
      max_coeff = TO_APNT(p2_coeffs[i]);
  }
  
  Sturm_seq<NT>  sts1 (p1);
  Sturm_seq<NT>  sts2 (p2);

  // Determine the number of iterations.
  const int      max_deg = (p1_deg > p2_deg) ? p1_deg : p2_deg;
  const int      max_iters = max_deg *
                             static_cast<int>(log(max_coeff) / log(2) + 1);
  NT             m1, m2;
  const NT       _two = NT(2);

  for (i = 0; i < max_iters; i++)
  {
    // Bisect [l1, r1] and select either [l1, m1] or [m1, r1].
    m1 = LEDA_NAMESPACE_NAME::to_bigfloat((r1 + l1) / _two);

    if (sts1.sign_changes_at_x(l1) - sts1.sign_changes_at_x(m1) == 1)
      r1 = m1;
    else
      l1 = m1;

    // Bisect [l2, r2] and select either [l2, m2] or [m2, r2].
    m2 = LEDA_NAMESPACE_NAME::to_bigfloat((r2 + l2) / _two);

    if (sts2.sign_changes_at_x(l2) - sts2.sign_changes_at_x(m2) == 1)
      r2 = m2;
    else
      l2 = m2;

    // If [l1, r1] and [l2, r2] do not overlap, return the comparison result.
    if (r1 < l2)
      return (SMALLER);
    else if (l1 > r2)
      return (LARGER);
  }

  // According to Loos, if we reached here, the two roots are equal.
  return (EQUAL);
}

CGAL_END_NAMESPACE

#endif
