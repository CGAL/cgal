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
// file          : include/CGAL/Arrangement_2/Polynom.h
// package       : Arrangement (2.62)
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_POLYNOM_H
#define CGAL_POLYNOM_H

#include <vector>
#include <iostream>

template <class NT>
class Polynom
{
 private:

  typename std::vector<NT> coeffs;          // The vector of coefficients.
  int                      degree;          // The degree of the polynomial.

 public:

  // Constructors.
  Polynom () :
    coeffs(),
    degree(-1)
  {}

  // Destructor.
  virtual ~Polynom ()
  {}

  // Get the degree of the polynomial (-1 means the polynomial equals 0).
  int deg () const
  {
    return (degree);
  }

  // Get the i'th coefficient of the polynomial.
  const NT& get (const int& i) const
  {
    return (coeffs[i]);
  }

  // Set the i'th coefficient of the polynomial.
  void set (const int& i, const NT& c)
  {
    // Resize the coefficients vector, if needed.
    if (i > degree)
    {
      coeffs.resize(i+1);
      for (int j = degree + 1; j < i; j++)
	coeffs[j] = NT(0);
      degree = i;
    }

    // Set the coefficient.
    coeffs[i] = c;

    // If needed, update the degree.
    _reduce_degree();

    return;
  }

  // Evaluate p(x0) (where p(x) is our polynomial).
  NT eval_at (const NT& x0) const
  {
    // Check the case of a zero polynomial.
    if (degree < 0)
      return (NT(0));

    // Use Horner's method:
    //  p(x) = a[n]*x^n + a[n-1]*x^(n-1) + ... + a[1]*x + a[0] =
    //       = (( ... (a[n]*x + a[n-1])*x + ... )*x + a[1])*x + a[0]
    NT     result = coeffs[degree];

    for (int i = degree - 1; i >= 0; i--)
      result = result*x0 + coeffs[i];

    return (result);
  }

  // Add two polynomials.
  Polynom<NT> operator+ (const Polynom<NT>& p) const
  {
    // Find the polynomial with the higher degree.
    const Polynom<NT>  *hdeg_P, *ldeg_P;

    if (degree > p.degree)
    {
      hdeg_P = this;
      ldeg_P = &p;
    }
    else
    {
      hdeg_P = &p;
      ldeg_P = this;
    }

    // Start with a copy of the polynomial of the higher degree.
    Polynom<NT>        res (*hdeg_P);
    
    // Add the coefficients of the other polynomial.
    for (int j = 0; j <= ldeg_P->degree; j++)
    {
      res.coeffs[j] += ldeg_P->coeffs[j];
    }

    // If needed, reduce the degree of the resulting polynomial.
    res._reduce_degree();

    // Return the result.
    return (res);
  }

  // Subtract two polynomials.
  Polynom<NT> operator- (const Polynom<NT>& p) const
  {
    // Find the polynomial with the higher degree.
    const Polynom<NT>  *hdeg_P, *ldeg_P;
    int                sign;

    if (degree > p.degree)
    {
      hdeg_P = this;
      ldeg_P = &p;
      sign = 1;
    }
    else
    {
      hdeg_P = &p;
      ldeg_P = this;
      sign = -1;
    }

    // Start with a copy of the polynomial of the higher degree.
    Polynom<NT>        res (hdeg_P->degree);
    int                j;

    for (j = 0; j <= hdeg_P->degree; j++)
    {
      if (sign == 1)
	res.coeffs[j] = hdeg_P->coeffs[j];
      else
	res.coeffs[j] = -(hdeg_P->coeffs[j]);
    }

    // Subtract the coefficients of the other polynomial.
    for (int j = 0; j <= ldeg_P->degree; j++)
    {
      if (sign == 1)
	res.coeffs[j] -= ldeg_P->coeffs[j];
      else
	res.coeffs[j] += ldeg_P->coeffs[j];
    }

    // If needed, reduce the degree of the resulting polynomial.
    res._reduce_degree();

    // Return the result.
    return (res);
  }

  // Multiply two polynomials.
  Polynom<NT> operator* (const Polynom<NT>& p) const
  {
    // In case one of the polynomial is zero, return a zero polynomial.
    if (degree < 0 || p.degree < 0)
    {
      Polynom<NT>   zero_poly;
      return (zero_poly);
    }

    // Allocate a polynomial whose degrees is the sum of the two degrees
    // and initialize it with zeroes.
    int         sum_deg = degree + p.degree;
    Polynom<NT> res (sum_deg);
    int         i, j;
    const NT    _zero (0);

    for (i = 0; i <= sum_deg; i++)
      res.coeffs[i] = _zero;
    
    // Perform the convolution.
    for (i = 0; i <= degree; i++)
      for (j = 0; j <= p.degree; j++)
	res.coeffs[i + j] += (coeffs[i] * p.coeffs[j]);

    // Return the result (no need for degree reduction).
    return (res);
  }

  // Divide two polynomials.
  Polynom<NT> operator/ (const Polynom<NT>& p) const
  {
    Polynom<NT> q, r;

    // Find q(x), r(x) such that:
    //  (*this)(x) = q(x)*p(x) + r(x)
    _divide (p, q, r);

    // Return the quontient polynomial.
    return (q);
  }

   // Return the residue from the division of the two polynomials.
  Polynom<NT> operator% (const Polynom<NT>& p) const
  {
    Polynom<NT> q, r;

    // Find q(x), r(x) such that:
    //  (*this)(x) = q(x)*p(x) + r(x)
    _divide (p, q, r);

    // Return the residue polynomial.
    return (r);
  }

  // Unary minus.
  Polynom<NT> operator- () const
  {
    // Negate all coefficients.
    Polynom<NT> res (degree);

    for (int i = 0; i <= degree; i++)
      res.coeffs[i] = -(coeffs[i]);

    return (res);
  }

  // Return the derivative of the polynomial.
  Polynom<NT> derive () const
  {
    // In case of a constant (or a zero) polynomial, the derivative is zero.
    if (degree <= 0)
    {
      Polynom<NT>   zero_poly;
      return (zero_poly);
    }

    // If out polynomial is:
    //  p(x) = a[n]*x^n + a[n-1]*x^(n-1) + ... + a[1]*x + a[0]
    // Then its derivative is:
    //  p(x) = n*a[n]*x^(n-1) + (n-1)*a[n-1]*x^(n-2) + ... + 2*a[2]*x + a[1]
    Polynom<NT> res (degree - 1);

    for (int i = 1; i <= degree; i++)
      res.coeffs[i - 1] = NT(i)*coeffs[i];

    // Return the result (no need for degree reduction).
    return (res);
  }

 protected:

  // Protected constructor: set the initial size of the coefficients vector.
  Polynom (const int& _degree) :
    coeffs(_degree + 1),
    degree(_degree)
  {}

  // Reduce the degree of the polynomial (omit zero coefficients).
  void _reduce_degree ()
  {
    int    old_deg = degree;

    // Omit any leading zero coefficients.
    while (degree >= 0 && coeffs[degree] == NT(0))
      degree--;

    // If the degree has been reduced, resize the coefficients vector.
    if (degree != old_deg)
      coeffs.resize(degree + 1);

    return;
  }

  // Calculate the quontient and residue polynomials of the division,
  // i.e. find q(x) and r(x) such that:
  //  (*this)(x) = q(x)*p(x) + r(x)
  void _divide (const Polynom<NT>& p,
		Polynom<NT>& q, Polynom<NT>& r) const
  {
    // In case p(x) has a higher degree, then q(x) is zero and r(x) is (*this).
    if (p.degree > degree)
    {
      Polynom<NT>   zero_poly;
      q = zero_poly;
      r = *this;
      return;
    }

    // Otherwise, q(x)'s degree is the difference between the two degrees, and
    // the final degree of r(x) should be less than p(x)'s degree.
    const NT    _zero (0);
    int         i;

    q.degree = degree - p.degree;
    q.coeffs.resize (q.degree + 1);

    for (i = 0; i <= q.degree; i++)
      q.coeffs[i] = _zero;

    // Perform "long division".
    Polynom<NT> work1(*this), work2;
    Polynom<NT> *work1_P = &work1, *work2_P = &work2, *swap_P;
    int         deg_diff;
    NT          factor;
    int         r_deg = degree;
    
    while (r_deg >= p.degree)
    {
      deg_diff = r_deg - p.degree;
      factor = work1_P->coeffs[r_deg] / p.coeffs[p.degree];

      q.coeffs[deg_diff] = factor;

      work2_P->coeffs.resize(r_deg);
      for (i = p.degree - 1; i >= 0; i--)
	work2_P->coeffs[i + deg_diff] = work1_P->coeffs[i + deg_diff] - 
	                                (factor * p.coeffs[i]);
      for (i = 0; i < deg_diff; i++)
	work2_P->coeffs[i] = work1_P->coeffs[i];
      
      do
      {
	r_deg--;
      } while (r_deg >= 0 && work2_P->coeffs[r_deg] == _zero);
      
      swap_P = work1_P;
      work1_P = work2_P;
      work2_P = swap_P;
    }

    // Assign the residue polynomial.
    if (r_deg < 0)
    {
      Polynom<NT>   zero_residue;
      r = zero_residue;
    }
    else
    {
      r.degree = r_deg;
      r.coeffs.resize (r.degree + 1);

      for (i = 0; i <= r.degree; i++)
	r.coeffs[i] = work1_P->coeffs[i];
    }

    return;
  }

};

// Print the polynomial.
template <class NT>
std::ostream& operator<< (std::ostream& os, const Polynom<NT>& p)
{
  const NT _zero(0);
  int      degree = p.deg();

  if (degree < 0)
  {
    os << '(' << _zero << ')';
    return (os);
  }

  os << '(' << p.get(degree) << ')';
  if (degree == 1)
    os << "*x";
  else if (degree > 1)
    os << "*x^" << degree;

  for (int i = degree - 1; i >= 0; i--)
  {
    const NT& c = p.get(i);

    if (c == 0)
      continue;
      os << " + (" << c << ')';
    
    if (i == 1)
      os << "*x";
    else if (i > 1)
      os << "*x^" << i;
  }

  return (os);
}

// Calculate the GCD of the two polynomials.
template <class NT>
Polynom<NT> polynom_gcd (const Polynom<NT>& f, const Polynom<NT>& g)
{
  // Let p0(x) be the polynomial with the higher degree, and p1(x)
  // be the polynomial with a lower degree.
  Polynom<NT> p0, p1;
  
  if (f.deg() > g.deg())
  {
    p0 = f;
    p1 = g; 
  }
  else
  {
    p0 = g;
    p1 = f;
  }
    
  // Use Euclid's algorithm.
  Polynom<NT>       r;

  do
  {
    r = p0 % p1;
      
    p0 = p1;
    p1 = r;
  } while (r.deg() > 0);

  if (r.deg() < 0)
  {
    // f(x), g(x) are not disjoint, and p0(x) is their GCD.
    return (p0);
  }
  else
  {
    // f(x), g(x) are disjoint, therefore their GCD is 1.
    Polynom<NT>  unit_gcd;
    unit_gcd.set (0, NT(1));

    return (unit_gcd);
  }
}

#endif
