// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Interval_arithmetic.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================


// This file contains the description of the two classes:
// - CGAL_Interval_nt_advanced  (do the FPU rounding mode changes yourself)
// - CGAL_Interval_nt		("plug-in" version, derived from the other one)
// The second one is slower.

// TODO:  operators for double * Interval -> Interval,
//                      Interval * double -> Interval
//        maybe a function (double * double) -> Interval.

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

#include <iostream.h>
#include <CGAL/assertions.h>
#include <CGAL/double.h>	// For CGAL_is_valid() and CGAL_is_finite().
#include <CGAL/Interval_arithmetic/_FPU.h>	// FPU rounding mode functions.


class CGAL_Interval_nt_advanced
{
  friend CGAL_Interval_nt_advanced sqrt(const CGAL_Interval_nt_advanced&);
  friend double CGAL_to_double(const CGAL_Interval_nt_advanced&);
  friend ostream& operator<<(ostream& os, const CGAL_Interval_nt_advanced&);

public:

  // The constructors.
  CGAL_Interval_nt_advanced() {}

  CGAL_Interval_nt_advanced(const double d)
  {
    inf = -d;
    sup = d;
  }

  CGAL_Interval_nt_advanced(const double i, const double s)
  {
#ifndef CGAL_NO_PRECONDITIONS
    CGAL_assertion(i<=s);
#endif
    inf = -i;
    sup = s;
  }

  CGAL_Interval_nt_advanced& operator=(const CGAL_Interval_nt_advanced& d)
  {
    inf = d.inf;
    sup = d.sup;
    return *this;
  }

  // The operators.
  CGAL_Interval_nt_advanced operator+(const CGAL_Interval_nt_advanced& d) const
  {
    return CGAL_Interval_nt_advanced(-(inf + d.inf), sup + d.sup);
  }

  CGAL_Interval_nt_advanced operator-(const CGAL_Interval_nt_advanced& d) const
  {
    return CGAL_Interval_nt_advanced(-(inf + d.sup), sup + d.inf);
  }

  CGAL_Interval_nt_advanced operator*(const CGAL_Interval_nt_advanced& d) const
  {
    if (inf<=0)					/* this>=0 */
    {
      if (d.inf<=0)				/* d>=0 */
	return CGAL_Interval_nt_advanced(-((-inf)*d.inf), sup*d.sup);
      else if (d.sup<=0)			/* d<=0 */
	return CGAL_Interval_nt_advanced(-(sup*d.inf), (-inf)*d.sup);
      else					/* 0 \in d */
	return CGAL_Interval_nt_advanced(-(sup*d.inf), sup*d.sup);
    }
    else if (sup<=0)				/* this<=0 */
    {
      if (d.inf<=0)				/* d>=0 */
	return CGAL_Interval_nt_advanced(-(inf*d.sup), sup*(-d.inf));
      else if (d.sup<=0)			/* d<=0 */
	return CGAL_Interval_nt_advanced(-((-sup)*d.sup), inf*d.inf);
      else					/* 0 \in d */
	return CGAL_Interval_nt_advanced(-(inf*d.sup), inf*d.inf);
    }
    else					/* 0 \in [inf;sup] */
    {
      if (d.inf<=0)				/* d>=0 */
	return CGAL_Interval_nt_advanced(-(inf*d.sup), sup*d.sup);
      else if (d.sup<=0)			/* d<=0 */
	return CGAL_Interval_nt_advanced(-(sup*d.inf), inf*d.inf);
      else					/* 0 \in d */
      {
        double tmp1, tmp2, tmp3, tmp4;
        tmp1 = inf*d.sup;
        tmp2 = sup*d.inf;
        tmp3 = inf*d.inf;
        tmp4 = sup*d.sup;
	return CGAL_Interval_nt_advanced(-((tmp1>tmp2) ? tmp1 : tmp2),
					   (tmp3>tmp4) ? tmp3 : tmp4);
      };
    };
  }

  CGAL_Interval_nt_advanced operator/(const CGAL_Interval_nt_advanced& d) const
  {
    if (d.inf<0.0)				/* d>0 */
    {
      if (inf<=0.0)				/* this>=0 */
	return CGAL_Interval_nt_advanced(-(inf/d.sup), sup/(-d.inf));
      else if (sup<=0.0)			/* this<=0 */
	return CGAL_Interval_nt_advanced(-(inf/(-d.inf)), sup/d.sup);
      else					/* 0 \in this */
	return CGAL_Interval_nt_advanced(-(inf/(-d.inf)), sup/(-d.inf));
    }
    else if (d.sup<0.0)				/* d<0 */
    {
      if (inf<=0.0)				/* this>=0 */
	return CGAL_Interval_nt_advanced(-(sup/(-d.sup)), inf/d.inf);
      else if (sup<=0.0)			/* b<=0 */
	return CGAL_Interval_nt_advanced(-(sup/d.inf), inf/(-d.sup));
      else					/* 0 \in this */
	return CGAL_Interval_nt_advanced(-(sup/(-d.sup)), inf/(-d.sup));
    }
    else					/* 0 \in [d.inf;d.sup] */
      return CGAL_Interval_nt_advanced(-(HUGE_VAL), HUGE_VAL);
  }

  CGAL_Interval_nt_advanced& operator+=(const CGAL_Interval_nt_advanced& d)
  {
    return *this = *this + d;
  }

  CGAL_Interval_nt_advanced& operator-=(const CGAL_Interval_nt_advanced& d)
  {
      // Check if this compact "one line" notation is ok for speed.
    return *this = *this - d;
    // inf += d.sup;
    // sup += d.inf;
    // return *this;
  }

  CGAL_Interval_nt_advanced & operator*=(const CGAL_Interval_nt_advanced& d)
  {
    return *this = *this * d;
    // return *this;
  }

  CGAL_Interval_nt_advanced & operator/=(const CGAL_Interval_nt_advanced& d)
  {
    return *this = *this / d;
    // return *this;
  }

  CGAL_Interval_nt_advanced operator-() const
  {
    return CGAL_Interval_nt_advanced(-(sup), inf);
  }

  bool is_valid() const
  {
    return CGAL_is_valid(inf) && CGAL_is_valid(sup) && (sup >= -inf);
  }

  bool is_finite() const
  {
    return CGAL_is_finite(inf) && CGAL_is_finite(sup);
  }

  bool operator==(const CGAL_Interval_nt_advanced& d) const
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
    CGAL_warning_msg(!overlap(d), " Comparison between overlapping intervals");
#endif
    return overlap(d);
  }

  bool operator!=(const CGAL_Interval_nt_advanced& d) const
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
    CGAL_warning_msg(!overlap(d), " Comparison between overlapping intervals");
#endif
    return !overlap(d);
  }

  bool operator<(const CGAL_Interval_nt_advanced& d) const
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
    CGAL_warning_msg(!overlap(d), " Comparison between overlapping intervals");
#endif
    return (sup < -d.inf);
  }

  bool operator>(const CGAL_Interval_nt_advanced& d) const
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
    CGAL_warning_msg(!overlap(d), " Comparison between overlapping intervals");
#endif
    return (d.sup < -inf);
  }

  bool operator<=(const CGAL_Interval_nt_advanced& d) const
  {
    return !(*this > d);
  }

  bool operator>=(const CGAL_Interval_nt_advanced& d) const
  {
    return !(*this < d);
  }

  double lower_bound() const { return -inf; }
  double upper_bound() const { return sup; }

private:
  bool overlap(const CGAL_Interval_nt_advanced& d) const
  {
    return (sup >= -d.inf) && (d.sup >= -inf);
  }

protected:
        // "inf" stores the OPPOSITE of the lower bound.
        // "sup" stores the upper bound of the interval.
  double inf, sup;
};


class CGAL_Interval_nt : public CGAL_Interval_nt_advanced
{
  friend CGAL_Interval_nt sqrt(const CGAL_Interval_nt&);

public:

  // Constructors are identical.
  CGAL_Interval_nt()
      {}
  CGAL_Interval_nt(const double d)
      : CGAL_Interval_nt_advanced(d) {}
  CGAL_Interval_nt(const double a, const double b)
      : CGAL_Interval_nt_advanced(a,b) {}
  CGAL_Interval_nt(const CGAL_Interval_nt_advanced &d)
      : CGAL_Interval_nt_advanced(d) {}

  // The member functions that have to be protected against rounding mode.
  CGAL_Interval_nt operator+(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator-(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator*(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator/(const CGAL_Interval_nt& d) const ;
  // The following do not seem to need a redefinition.
#if 0
  CGAL_Interval_nt& operator+=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator-=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator*=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator/=(const CGAL_Interval_nt& d) ;
#endif
};


inline CGAL_Interval_nt CGAL_Interval_nt::operator+(const CGAL_Interval_nt& d)
    const // return tmp;
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator+(d) );
  // tmp = CGAL_Interval_nt_advanced::operator+(d);
  // CGAL_Interval_nt tmp (-(inf + d.inf), sup + d.sup);
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator-(const CGAL_Interval_nt& d)
    const
{
  CGAL_FPU_set_rounding_to_infinity();
  // CGAL_Interval_nt tmp (-(inf + d.sup), sup + d.inf);
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator-(d) );
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator*(const CGAL_Interval_nt& d)
    const
{
  CGAL_FPU_set_rounding_to_infinity();
  // CGAL_Interval_nt tmp = ((CGAL_Interval_nt_advanced) *this)
  //	* (CGAL_Interval_nt_advanced) d;
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator*(d) );
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator/(const CGAL_Interval_nt& d)
    const
{
  CGAL_FPU_set_rounding_to_infinity();
  // CGAL_Interval_nt tmp = ((CGAL_Interval_nt_advanced) *this)
  //	/ (CGAL_Interval_nt_advanced) d;
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator/(d) );
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

#if 0
inline CGAL_Interval_nt& CGAL_Interval_nt::operator+=(const CGAL_Interval_nt& d)
{
    //  Stress/Bench test this approach.
    //  Right now, the test-suite doesn't cover this case.
  return *this = *this + d;

    /*
  CGAL_FPU_set_rounding_to_infinity();
  // inf += d.inf; sup += d.sup;
  CGAL_Interval_nt_advanced::operator+=(d);
  CGAL_FPU_set_rounding_to_nearest();
  return *this;
  */
}

inline CGAL_Interval_nt& CGAL_Interval_nt::operator-=(const CGAL_Interval_nt& d)
{
  CGAL_FPU_set_rounding_to_infinity();
  // inf += d.sup; sup += d.inf;
  CGAL_Interval_nt_advanced::operator-=(d);
  CGAL_FPU_set_rounding_to_nearest();
  return *this;
}

inline CGAL_Interval_nt& CGAL_Interval_nt::operator*=(const CGAL_Interval_nt& d)
{
  CGAL_FPU_set_rounding_to_infinity();
  // *this = ((CGAL_Interval_nt_advanced) *this)
  //	* (CGAL_Interval_nt_advanced) d;
  CGAL_Interval_nt_advanced::operator*=(d);
  CGAL_FPU_set_rounding_to_nearest();
  return *this;
}

inline CGAL_Interval_nt& CGAL_Interval_nt::operator/=(const CGAL_Interval_nt& d)
{
  CGAL_FPU_set_rounding_to_infinity();
  // *this = ((CGAL_Interval_nt_advanced) *this)
  //	/ (CGAL_Interval_nt_advanced) d;
  CGAL_Interval_nt_advanced::operator/=(d);
  CGAL_FPU_set_rounding_to_nearest();
  return *this;
}
#endif

inline CGAL_Interval_nt_advanced sqrt(const CGAL_Interval_nt_advanced& d)
{
  CGAL_FPU_set_rounding_to_minus_infinity();
  CGAL_Interval_nt_advanced tmp;
  tmp.inf = - sqrt(-d.inf);
  CGAL_FPU_set_rounding_to_infinity();
  tmp.sup = sqrt(d.sup);
  return tmp;
}

inline CGAL_Interval_nt sqrt(const CGAL_Interval_nt& d) // return tmp;
{
  CGAL_Interval_nt tmp;
  tmp = sqrt( (CGAL_Interval_nt_advanced) d);
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

// This one need a version for ..._advanced ?
// Is this function only needed by GPC...
inline CGAL_Interval_nt CGAL_to_interval_nt(const double &d)
{
    return (CGAL_Interval_nt) d;
}

inline double CGAL_to_double(const CGAL_Interval_nt_advanced& d)
{
    return (d.sup-d.inf)*.5;
}

ostream& operator<<(ostream& os, const CGAL_Interval_nt_advanced& d)
{
  return os << "[" << -d.inf << ";" << d.sup << "]";
}


// Finally we source the cast functions from other NTs, when necessary.

#ifdef CGAL_GMPZ_H
#include <CGAL/Interval_arithmetic/IA_Gmpz.h>
#endif

#ifdef CGAL_BIGFLOAT_H
#include <CGAL/Interval_arithmetic/IA_leda_bigfloat.h>
#endif

#ifdef CGAL_INTEGER_H
#include <CGAL/Interval_arithmetic/IA_leda_integer.h>
#endif

#ifdef CGAL_REAL_H
#include <CGAL/Interval_arithmetic/IA_leda_real.h>
#endif

#ifdef CGAL_RATIONAL_H
#include <CGAL/Interval_arithmetic/IA_leda_rational.h>
#endif

#endif /* CGAL_INTERVAL_ARITHMETIC_H */
