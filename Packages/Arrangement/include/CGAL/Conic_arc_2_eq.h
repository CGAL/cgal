#ifndef CGAL_CONIC_ARC_2_EQ_H
#define CGAL_CONIC_ARC_2_EQ_H

//#include <CGAL/number_type_tags.h>
//#include <CGAL/double.h>
//#include <CGAL/double.h>
#include <CGAL/basic.h>

#include <CGAL/leda_bigfloat.h>

CGAL_BEGIN_NAMESPACE

/*
typedef double APNT;
#define APNT_EPS    0.00000001
#define COMP_EPS    0.000001
#define TO_APNT(x)  CGAL::to_double(x)
#define APNT_ABS(x) fabs(x)
*/

typedef leda_bigfloat  APNT;
#define APNT_BITLEN    512
#define ZERO_BITLEN    -25
#define APNT_ZERO(x)   (x == 0 || ilog2(x) < ZERO_BITLEN)
#define APNT_CZERO(x)  (x == 0 || ilog2(x) < ZERO_BITLEN/2)
#define TO_APNT(x)     to_bigfloat(x)
#define APNT_ABS(x)    CGAL_NTS abs(x)


// ----------------------------------------------------------------------------
// Compare two approximate values.
//
template <class APNT>
Comparison_result eps_compare (const APNT& val1, const APNT& val2)
{
  APNT denom = APNT_ABS(val1 - val2);
  APNT numer = APNT_ABS(val1) + APNT_ABS(val2);

  if (APNT_CZERO(val1) && APNT_CZERO(val2))
  //  if (numer < APNT_EPS)
  {
    // The two number are very close to zero:
    if (APNT_ZERO(val1 - val2))
    //if (denom < numer)
      return (EQUAL);
    else
      return ((val1 < val2) ? SMALLER : LARGER);
  }
  else
  {
    //                  | x - y |
    //  x = y    <==>  ----------- < COMP_EPS
    //                  |x| + |y|
    if (APNT_ZERO(denom / numer))
    //if (denom / numer < COMP_EPS)
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
		    NT& sqrt_re, NT& sqrt_im)
{
  static const NT    _Zero = 0;
  static const NT    _Two = 2;

  // Check the special case where imag=0:
  if (imag == _Zero)
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
			 NT *roots)
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

  if (R_square == _Zero)
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
		     D_re, D_im);
      _complex_sqrt (E_square_re, E_square_im,
		     E_re, E_im);

      // Use only the combinations that yield real-valued numbers.
      if (R_im == -D_im)
      {
	roots[n_roots] = -_One_quarter*a3 + _One_half*D_re;
	n_roots++;
      }
      if (R_im == D_im)
      {
	roots[n_roots] = -_One_quarter*a3 - _One_half*D_re;
	n_roots++;
      }

      if (R_im == E_im)
      {
	roots[n_roots] = -_One_quarter*a3 + _One_half*E_re;
	n_roots++;
      }
      if (R_im == -E_im)
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
			NT* roots)
{
  static const NT _zero = 0;
  static const NT _two = 2;
  static const NT _four = 4;

  if (a == _zero)
  {
    if (b != _zero)
    {
      roots[0] = -c/b;
      return (1);
    }
    else
      return (0);
  }

  const NT        disc = b*b - _four*a*c;

  if (disc < _zero)
  {
    return (0);
  }
  else if (disc == _zero)
  {
    roots[0] = roots[1] = -b/(_two*a);
    return (1);
  }
  else
  {
    const NT      sqrt_disc = CGAL::sqrt(disc);

    roots[0] = (sqrt_disc - b)/(_two*a);
    roots[1] = -(sqrt_disc + b)/(_two*a);
    return (2);
  }
}

template <class NT>
int _solve_quadratic_eq (const NT& a, const NT& b, const NT& c,
			 NT* roots)
{
  static const NT _zero = 0;
  static const NT _two = 2;
  static const NT _four = 4;

  if (a == _zero)
  {
    if (b != _zero)
    {
      roots[0] = -c/b;
      return (1);
    }
    else
      return (0);
  }

  const NT        disc = b*b - _four*a*c;

  if (disc < _zero)
  {
    return (0);
  }
  else if (disc == _zero)
  {
    roots[0] = roots[1] = -b/(_two*a);
    return (1);
  }
  else
  {
    const NT      sqrt_disc = sqrt_d(disc, APNT_BITLEN, 2);

    roots[0] = (sqrt_disc - b)/(_two*a);
    roots[1] = -(sqrt_disc + b)/(_two*a);
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
			    NT* roots, int& n_approx)
{
  n_approx = 0;

  // In case we have an equation of a lower degree:
  static const NT _zero = 0;

  if (_a == _zero)
  {
    // We have to solve _b*x^2 + _c*x + _d = 0.
    return (solve_quadratic_eq<NT> (_b, _c, _d,
				    roots));
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
    // We can obtainf this root using the fact that -b is the sum of p(x)'s 
    // roots.
    if (B == _zero)
    {
      roots[0] = -b / _three;
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
      roots[1] = -(b + 2*x0);
      return (2);
    }
  }
      
  // If we reached here, we can be sure that there are no multiple roots.
  // We continue by trying to find a quadratic factor (x^2 - r*x - q) that
  // divides p(x).
  const APNT        b_app = TO_APNT(b);
  const APNT        c_app = TO_APNT(c);
  const APNT        d_app = TO_APNT(d);
  APNT              q = 1;
  APNT              r = 2;
  APNT              A0, B0;
  APNT              A1, B1;
  APNT              mat[2][2];
  APNT              det;
  APNT              r_tag, q_tag;
  bool              another_iter = true;
  int               iter = 0;

  while (another_iter && iter < APNT_BITLEN)
  {
    iter++;

    // Let us denote:
    //  p(x) = p1(x) * (x^2 - r*x - q) + (A0*x + B0)
    //  p1(x) = 0 * (x^2 - r*x - q) + (A1*x + B1)
    //
    // Where:
    //  p1(x) = x + psi
    //
    // Then we can obtain:
    A0 = c_app + q + r*(b_app + r);
    B0 = d_app + q*(b_app + r);
    A1 = 1;
    B1 = b_app + r;

    // Now exchange r,q with the updated approximations r',q' given by:
    //
    //  +-  -+   +- -+   +-          -+ -1 +-  -+
    //  | r' | = | r | - | A1*r+B1 A1 |    | A0 |
    //  | q' |   | q |   |  A1*q   B1 |    | B0 |
    //  +-  -+   +- -+   +-          -+    +-  -+
    //
    mat[0][0] = A1*r + B1;
    mat[0][1] = A1;
    mat[1][0] = A1*q;
    mat[1][1] = B1;

    det = mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];
    if (APNT_ZERO(det))
    //    if (APNT_ABS(det) < APNT_EPS)
      break;

    r_tag = r - (mat[1][1]*A0 - mat[0][1]*B0) / det;
    q_tag = q - (-mat[1][0]*A0 + mat[0][0]*B0) / det;

    if (eps_compare<APNT> (r, r_tag) == EQUAL &&
	eps_compare<APNT> (q, q_tag) == EQUAL)
    {
      another_iter = false;
    }

    r = r_tag;
    q = q_tag;
  }

  // Now we have the approximation:
  //  p(x) = (x^2 + psi) * (x^2 - r*x - q)
  //
  // Try to solve the two quadratic factor:
  APNT       app_roots[2];

  n_approx = _solve_quadratic_eq<APNT> (APNT(1), -r, -q,
					app_roots + n_approx);
  
  // Copy the approximations.
  for (int i = 0; i < n_approx; i++)
    roots[i] = NT(app_roots[i]);

  // Add the extra root.
  roots[n_approx] = -(b_app + r);
  n_approx++;

  return (n_approx);
}

/*
template <class APNT>
APNT refine_solution (const APNT& a, const APNT& b, const APNT& c, 
		      const APNT& d, const APNT& e,
		      const APNT& _x)
{
  APNT    x = _x;
  APNT    val, deriv;
  
  for (int i = 0; i < 10; i++)
  {
    val = (((a*x + b)*x + c)*x + d)*x + e;

    if (APNT_ZERO(val))
      return (x);

    deriv = ((4*a*x + 3*b)*x + 2*c)*x + d;

    x -= val/deriv;
  }

  return (x);
}
*/

// ----------------------------------------------------------------------------
// Solve a quartic equation: _a*x^4 + _b*x^3 + _c*x^2 + _d*x + _e = 0.
// The function returns the number of distinct solutions.
// The roots area must be at least of size 4.
// 
template <class NT>
int solve_quartic_eq (const NT& _a, const NT& _b, const NT& _c, 
		      const NT& _d, const NT& _e,
		      NT* roots, int& n_approx)
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
	  return (1);
	}
	else
	{
	  // Add the solution 0 (with multiplicity 3), and add the solution
	  // to _a*x + _b = 0.
	  if (_a == _zero)
	    return (0);

	  roots[0] = _zero;
	  roots[1] = -(_b/_a);
	  return (2);
	}
      }
      else
      {
	// Add the solution 0 (with multiplicity 2),
	// and solve _a*x^2 + _b*x + _c*x = 0.
	roots[0] = _zero;

	return (1 + solve_quadratic_eq<NT> (_a, _b, _c,
					    roots + 1));
      }
    }
    else
    {
      // Add the solution 0, and solve _a*x^3 + _b*x^2 + _c*x + _d = 0.
      roots[0] = _zero;

      return (1 + _solve_cubic_eq<NT> (_a, _b, _c, _d,
				       roots + 1, n_approx));
    }
  }

  // In case we have an equation of a lower degree:
  if (_a == _zero)
  {
    if (_b == _zero)
    {
      // We have to solve _c*x^2 + _d*x + _e = 0.
      return (solve_quadratic_eq<NT> (_c, _d, _e,
				      roots));
    }
    else
    {
      // We have to solve _b*x^3 + _c*x^2 + _d*x + _e = 0.
      return (_solve_cubic_eq<NT> (_b, _c, _d, _e,
				   roots, n_approx));
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
    int  n_sq;
    int  n = 0;

    n_sq = solve_quadratic_eq<NT> (_a, _c, _e,
				   sq_roots);

    // Convert to roots of the original equation: only positive roots for y
    // are relevant of course.
    for (int j = 0; j < n_sq; j++)
    {
      if (sq_roots[j] > _zero)
      {
	roots[n] = CGAL::sqrt(sq_roots[j]);
	roots[n+1] = -roots[n];
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
	roots[1] = -(b + 3*x0);
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
					   roots);

	if (m == 1)
	{
	  // Add the second root.
	  roots[1] = -(b + 3*roots[0]);
	  m++;
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

	// The other two solutions are obtained from solving:
	// p(x) / (x - x0)^2 = x^2 + (b + 2*x0)*x + (c + 2*b*x0 + 3*x0^2)
	return (1 + 
		solve_quadratic_eq<NT> (NT(1), b + 2*x0, c + x0*(2*b + 3*x0),
					roots + 1));		
      }
    }
  }
  
  // If we reached here, we can be sure that there are no multiple roots.
  // We continue by trying to find a quadratic factor (x^2 - r*x - q) that
  // divides p(x).
  const APNT        b_app = TO_APNT(b);
  const APNT        c_app = TO_APNT(c);
  const APNT        d_app = TO_APNT(d);
  const APNT        e_app = TO_APNT(e);
  APNT              q = 1;
  APNT              r = 2;
  APNT              phi, psi;
  APNT              A0, B0;
  APNT              A1, B1;
  APNT              mat[2][2];
  APNT              det;
  APNT              r_tag, q_tag;
  bool              another_iter = true;
  int               iter = 0;

  double            droots[4];
  int               n_droots = _Ferrari_quartic_eq<double> (1,
							    CGAL::to_double(b),
							    CGAL::to_double(c),
							    CGAL::to_double(d),
							    CGAL::to_double(e),
							    droots);
  
  if (n_droots == 0)
    return (0);
  r = APNT(droots[0] + droots[1]);
  q = -APNT(droots[0] * droots[1]);
  
  while (another_iter && iter < APNT_BITLEN)
  {
    iter++;

    // Let us denote:
    //  p(x) = p1(x) * (x^2 - r*x - q) + (A0*x + B0)
    //  p1(x) = p2(x) * (x^2 - r*x - q) + (A1*x + B1)
    //
    // Where:
    //  p1(x) = x^2 + phi*x + psi
    //
    // Then we can obtain:
    A0 = d_app + b_app*q + r*(2*q + c_app + r*(b_app + r));
    B0 = e_app + q*(c_app + r*(b_app + r) + q);
    A1 = b_app + 2*r;
    B1 = c_app + 2*q + r*(b_app + r);

    // Now exchange r,q with the updated approximations r',q' given by:
    //
    //  +-  -+   +- -+   +-          -+ -1 +-  -+
    //  | r' | = | r | - | A1*r+B1 A1 |    | A0 |
    //  | q' |   | q |   |  A1*q   B1 |    | B0 |
    //  +-  -+   +- -+   +-          -+    +-  -+
    //
    mat[0][0] = A1*r + B1;
    mat[0][1] = A1;
    mat[1][0] = A1*q;
    mat[1][1] = B1;

    det = mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];
    if (APNT_ZERO(det))
      //if (APNT_ABS(det) < APNT_EPS)
      break;

    r_tag = r - (mat[1][1]*A0 - mat[0][1]*B0) / det;
    q_tag = q - (-mat[1][0]*A0 + mat[0][0]*B0) / det;

    if (eps_compare<APNT> (r, r_tag) == EQUAL &&
	eps_compare<APNT> (q, q_tag) == EQUAL)
    {
      another_iter = false;
    }

    r = r_tag;
    q = q_tag;
  }

  // Now we have the approximation:
  //  p(x) = (x^2 + phi*x^2 + psi) * (x^2 - r*x - q)
  //
  // Try to solve the two quadratic factors:
  APNT       app_roots[4];

  phi = b_app + r;
  psi = c_app + q + r*(b_app + r);

  n_approx = _solve_quadratic_eq<APNT> (APNT(1), phi, psi,
					app_roots);
  n_approx += _solve_quadratic_eq<APNT> (APNT(1), -r, -q,
					 app_roots + n_approx);
  
  // Copy the approximations.
  for (int i = 0; i < n_approx; i++)
    roots[i] = NT(app_roots[i]);
  //roots[i] = NT(refine_solution(TO_APNT(_a), TO_APNT(_b), TO_APNT(_c),
  //				  TO_APNT(_d), TO_APNT(_e),
  //                              app_roots[i]));

  return (n_approx);
}

CGAL_END_NAMESPACE

#endif
