// Copyright (c) 1999  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_STURM_SEQ_H
#define CGAL_STURM_SEQ_H

#include <CGAL/Arrangement_2/Polynom.h>
#include <vector>
#include <list>

template <class NT>
class Sturm_seq
{
 private:

  std::list<Polynom<NT> > fs;
  Sturm_seq               *next_P;

  // Copy constructor and assignment operator - not supported.
  Sturm_seq (const Sturm_seq<NT>& );
  void operator= (const Sturm_seq<NT>& );

 public:

  // Construct the Strum sequence for the given polynomial.
  Sturm_seq (const Polynom<NT>& p) :
    fs(),
    next_P(NULL)
  {
    // If the polynomial is contant, just push it to the list and finish.
    if (p.deg() <= 0)
    {
      fs.push_back(p);
      return;
    }

    // We shall start with f[0](x) = p(x), f[1](x) = p'(x).
    Polynom<NT>  f0 (p);
    Polynom<NT>  f1 = p.derive();

    if (f1.deg() > 0)
    {
      // Calculate the greatest common divisor of p(x) and p'(x).
      Polynom<NT> gcd = polynom_gcd<NT> (f0, f1);
    
      if (gcd.deg() >= 0)
      {
	// Because two polynomials p(x) and p'(x) are not disjoint, then p(x) 
	// has at least one root with multiplicity >= 1:
	// First, construct a sequence for the GCD. 
	next_P = new Sturm_seq<NT> (gcd);

	// We know that q(x) = p(x)/GCD(p(x),p'(x)) has only simple roots,
	// so we start the sequence with f[0](x) = q(x), f[1](x) = q'(x).
	f0 = p / gcd;
	f1 = f0.derive();
      }
    }

    fs.push_back(f0);
    fs.push_back(f1);

    // Complete the rest of the sequence.
    const Polynom<NT> *fi_P = &f0;
    const Polynom<NT> *fi1_P = &f1;
    Polynom<NT>       fi2;

    while (fi1_P->deg() > 0)
    {
      // The next polynomial f[i+2] is defined by -(f[i](x) % f[i+1](x)).
      fs.push_back (-((*fi_P) % (*fi1_P)));
      
      fi_P = fi1_P;
      fi1_P = &(fs.back());
    }
  }

  // Destructor.
  ~Sturm_seq ()
  {
    if (next_P != NULL)
      delete next_P;
  }

  // Get the number of sign-changes in the Sturm sequence for a given x.
  int sign_changes_at_x (const NT& x) const
  {
    // Count the sign changes in the sequence:
    //  f[0](x), f[1](x), ... , f[s](x).
    typename std::list<Polynom<NT> >::const_iterator iter;
    NT                                 y;
    const NT                           _zero (0);
    int                                prev_sign = 0;
    int                                curr_sign;
    int                                result = 0;

    for (iter = fs.begin(); iter != fs.end(); ++iter)
    {
      // Determine the sign of f[i](x).
      y = (*iter).eval_at(x);

      if (y == _zero)
	curr_sign = 0;
      else
	curr_sign = (y > _zero) ? 1 : -1;
      

      // Count if there was a sign change.
      if (prev_sign != 0 && curr_sign != prev_sign)
	result++;

      prev_sign = curr_sign;
    }
    
    // Sum up the result with the number of sign changes in the next sequence.
    if (next_P != NULL)
      return (result + next_P->sign_changes_at_x (x));
    else
      return (result);
  }

  // Get the number of sign-changes at Infinity (if inf > 0) or at -Infinity
  // (if inf < 0).
  int sign_changes_at_infinity (const int& inf) const
  {
    // Count the sign changes in the sequence:
    //  f[0](-Inf), f[1](-Inf), ... , f[s](-Inf).
    typename std::list<Polynom<NT> >::const_iterator iter;
    int                                deg;
    NT                                 leading_coeff;
    const NT                           _zero(0);
    int                                prev_sign = 0;
    int                                curr_sign;
    int                                result = 0;

    for (iter = fs.begin(); iter != fs.end(); ++iter)
    {
      // Determine the sign of f[i](-Infinity): This is the sign of the leading
      // coefficient if the polynomial has an even dergee, otherwise it has an
      // opposite sign.
      deg = (*iter).deg();
      leading_coeff = (*iter).get(deg);

      curr_sign = (leading_coeff < _zero) ? -1 : 1;
      if (inf < 0 && deg % 2 == 1)
	curr_sign = -curr_sign;

      // Count if there was a sign change.
      if (prev_sign != 0 && curr_sign != prev_sign)
	result++;

      prev_sign = curr_sign;
    }
            
    // Sum up the result with the number of sign changes in the next sequence.
    if (next_P != NULL)
      return (result + next_P->sign_changes_at_infinity(inf));
    else
      return (result);
  }
};

#endif
