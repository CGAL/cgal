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

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

#include <iostream.h>
#include <CGAL/assertions.h>
#include <CGAL/double.h>	// For CGAL_is_valid() and CGAL_is_finite().
#include <CGAL/number_utils.h>  // For CGAL_{sign,compare,abs,min,max...}
#include <CGAL/Interval_arithmetic/_FPU.h>	// FPU rounding mode functions.

struct CGAL_Interval_nt_advanced
{
  struct unsafe_comparison{};  // Exception class.
protected:
  double inf, sup;	// "inf" stores the OPPOSITE of the lower bound.
			// "sup" stores the upper bound of the interval.
private:
  int overlap_action() const throw (unsafe_comparison)
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
    CGAL_warning_msg(false, " Comparison between overlapping intervals");
#endif
#ifndef CGAL_IA_DONT_THROW_EXCEPTION
    throw unsafe_comparison();
#endif
    return 0; // Return a value to prevent compiler warnings.
  }

public:
  friend CGAL_Interval_nt_advanced sqrt     (const CGAL_Interval_nt_advanced&);
  friend CGAL_Interval_nt_advanced CGAL_abs (const CGAL_Interval_nt_advanced&);
  friend CGAL_Interval_nt_advanced CGAL_min (const CGAL_Interval_nt_advanced&,
	  				     const CGAL_Interval_nt_advanced&);
  friend CGAL_Interval_nt_advanced CGAL_max (const CGAL_Interval_nt_advanced&,
					     const CGAL_Interval_nt_advanced&);
  friend CGAL_Interval_nt_advanced operator* (const double,
					      const CGAL_Interval_nt_advanced&);
  friend double CGAL_to_double (const CGAL_Interval_nt_advanced&);
  friend bool   CGAL_is_valid  (const CGAL_Interval_nt_advanced&);
  friend bool   CGAL_is_finite (const CGAL_Interval_nt_advanced&);
  friend CGAL_Sign CGAL_sign   (const CGAL_Interval_nt_advanced&);
  friend CGAL_Comparison_result CGAL_compare (const CGAL_Interval_nt_advanced&,
					      const CGAL_Interval_nt_advanced&);

  // The constructors.
  CGAL_Interval_nt_advanced() {}

  CGAL_Interval_nt_advanced(const double d)
      : inf(-d), sup(d) {}

  CGAL_Interval_nt_advanced(const double i, const double s)
      : inf(-i), sup(s)
#ifdef CGAL_NO_PRECONDITIONS
      {}
#else
      { CGAL_assertion(i<=s); }
#endif

#if 1
  // The copy constructors/assignment: useless.
  // The default ones are ok, but these are faster...
  CGAL_Interval_nt_advanced(const CGAL_Interval_nt_advanced& d)
      : inf(d.inf), sup(d.sup) {}

  CGAL_Interval_nt_advanced& operator=(const CGAL_Interval_nt_advanced& d)
  { inf = d.inf; sup = d.sup; return *this; }
#endif

  // The operators.
  CGAL_Interval_nt_advanced operator+(const CGAL_Interval_nt_advanced& d) const
  { return CGAL_Interval_nt_advanced(-(inf + d.inf), sup + d.sup); }
  // { return CGAL_Interval_nt_advanced (d) += *this; }

  CGAL_Interval_nt_advanced operator-(const CGAL_Interval_nt_advanced& d) const
  { return CGAL_Interval_nt_advanced(-(inf + d.sup), sup + d.inf); }

  // Those 2 ones could be made not inlined.
  CGAL_Interval_nt_advanced operator*(const CGAL_Interval_nt_advanced& d) const;
  CGAL_Interval_nt_advanced operator/(const CGAL_Interval_nt_advanced& d) const;

  CGAL_Interval_nt_advanced operator-() const
  { return CGAL_Interval_nt_advanced(-(sup), inf); }

  CGAL_Interval_nt_advanced& operator+=(const CGAL_Interval_nt_advanced& d);
  CGAL_Interval_nt_advanced& operator-=(const CGAL_Interval_nt_advanced& d);
  CGAL_Interval_nt_advanced& operator*=(const CGAL_Interval_nt_advanced& d);
  CGAL_Interval_nt_advanced& operator/=(const CGAL_Interval_nt_advanced& d);

  // For speed...
  CGAL_Interval_nt_advanced operator*(const double d) const;

  bool operator<(const CGAL_Interval_nt_advanced& d) const
  {
    if (sup  < -d.inf) return true;
    if (-inf >= d.sup) return false;
    return overlap_action();
  }

  bool operator>(const CGAL_Interval_nt_advanced& d) const
  { return (d < *this); }

  bool operator<=(const CGAL_Interval_nt_advanced& d) const
  {
    if (sup <= -d.inf) return true;
    if (-inf >  d.sup) return false;
    return overlap_action();
  }

  bool operator>=(const CGAL_Interval_nt_advanced& d) const
  { return (d <= *this); }

  bool operator==(const CGAL_Interval_nt_advanced& d) const
  {
    if ((-d.inf >  sup) || (d.sup  < -inf)) return false;
    if ((-d.inf == sup) && (d.sup == -inf)) return true;
    return overlap_action();
  }

  bool operator!=(const CGAL_Interval_nt_advanced& d) const
  { return !(d == *this); }

  double lower_bound() const { return -inf; }
  double upper_bound() const { return sup; }
};


inline CGAL_Interval_nt_advanced CGAL_Interval_nt_advanced::operator*
  (const CGAL_Interval_nt_advanced& d) const
{
  if (inf<=0)					/* this>=0 */
  {
    if (d.inf<=0)				/* d>=0 */
      return CGAL_Interval_nt_advanced(-((-inf)*d.inf), sup*d.sup);
    else if (d.sup<=0)				/* d<=0 */
      return CGAL_Interval_nt_advanced(-(sup*d.inf), (-inf)*d.sup);
    else					/* 0 \in d */
      return CGAL_Interval_nt_advanced(-(sup*d.inf), sup*d.sup);
  }
  else if (sup<=0)				/* this<=0 */
  {
    if (d.inf<=0)				/* d>=0 */
      return CGAL_Interval_nt_advanced(-(inf*d.sup), sup*(-d.inf));
    else if (d.sup<=0)				/* d<=0 */
      return CGAL_Interval_nt_advanced(-((-sup)*d.sup), inf*d.inf);
    else					/* 0 \in d */
      return CGAL_Interval_nt_advanced(-(inf*d.sup), inf*d.inf);
  }
  else						/* 0 \in [inf;sup] */
  {
    if (d.inf<=0)				/* d>=0 */
      return CGAL_Interval_nt_advanced(-(inf*d.sup), sup*d.sup);
    else if (d.sup<=0)				/* d<=0 */
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

inline CGAL_Interval_nt_advanced CGAL_Interval_nt_advanced::operator*
  (const double d) const
{
  if (d>=0)	return CGAL_Interval_nt_advanced(-(inf*d), sup*d);
  else		return CGAL_Interval_nt_advanced(-(sup*(-d)), inf*(-d));
}

inline CGAL_Interval_nt_advanced operator*
  (const double d, const CGAL_Interval_nt_advanced &t)
{ return t*d; }

inline CGAL_Interval_nt_advanced CGAL_Interval_nt_advanced::operator/
  (const CGAL_Interval_nt_advanced& d) const
{
  if (d.inf<0.0)				/* d>0 */
  {
    if (inf<=0.0)				/* this>=0 */
      return CGAL_Interval_nt_advanced(-(inf/d.sup), sup/(-d.inf));
    else if (sup<=0.0)				/* this<=0 */
      return CGAL_Interval_nt_advanced(-(inf/(-d.inf)), sup/d.sup);
    else					/* 0 \in this */
      return CGAL_Interval_nt_advanced(-(inf/(-d.inf)), sup/(-d.inf));
  }
  else if (d.sup<0.0)				/* d<0 */
  {
    if (inf<=0.0)				/* this>=0 */
      return CGAL_Interval_nt_advanced(-(sup/(-d.sup)), inf/d.inf);
    else if (sup<=0.0)				/* this<=0 */
      return CGAL_Interval_nt_advanced(-(sup/d.inf), inf/(-d.sup));
    else					/* 0 \in this */
      return CGAL_Interval_nt_advanced(-(sup/(-d.sup)), inf/(-d.sup));
  }
  else						/* 0 \in [d.inf;d.sup] */
    return CGAL_Interval_nt_advanced(-HUGE_VAL, HUGE_VAL);
	   // We could do slightly better -> [0;HUGE_VAL] when d.sup==0,
	   // but is this worth ?
}

inline CGAL_Interval_nt_advanced& CGAL_Interval_nt_advanced::operator+=
  (const CGAL_Interval_nt_advanced& d)
{
  // A pity: this compact "one line" notation is not ok for speed.
  // return *this = *this + d;
  inf += d.inf;
  sup += d.sup;
  return *this;
}

inline CGAL_Interval_nt_advanced& CGAL_Interval_nt_advanced::operator-=
  (const CGAL_Interval_nt_advanced& d)
{ // return *this = *this - d;
  inf += d.sup;
  sup += d.inf;
  return *this;
}

inline CGAL_Interval_nt_advanced& CGAL_Interval_nt_advanced::operator*=
  (const CGAL_Interval_nt_advanced& d)
{ return *this = *this * d; }

inline CGAL_Interval_nt_advanced& CGAL_Interval_nt_advanced::operator/=
  (const CGAL_Interval_nt_advanced& d)
{ return *this = *this / d; }

inline CGAL_Interval_nt_advanced sqrt(const CGAL_Interval_nt_advanced& d)
{
  CGAL_FPU_set_rounding_to_minus_infinity();
  CGAL_Interval_nt_advanced tmp;
  tmp.inf = - sqrt(-d.inf);
  CGAL_FPU_set_rounding_to_infinity();
  tmp.sup = sqrt(d.sup);
  return tmp;
}

inline double CGAL_to_double (const CGAL_Interval_nt_advanced& d)
{ return (d.sup-d.inf)*.5; }

inline bool CGAL_is_valid  (const CGAL_Interval_nt_advanced& d)
{ return CGAL_is_valid(d.inf) && CGAL_is_valid(d.sup) && (d.sup >= -d.inf); }

inline bool CGAL_is_finite (const CGAL_Interval_nt_advanced& d)
{ return CGAL_is_finite(d.inf) && CGAL_is_finite(d.sup); }

inline CGAL_Sign CGAL_sign (const CGAL_Interval_nt_advanced& d)
{
    // Benchmark it, and compare with CGAL_compare(d,0).
  if (d.inf < 0) return CGAL_POSITIVE;
  if (d.sup < 0) return CGAL_NEGATIVE;
  if ( (d.inf == 0) && (d.sup == 0) ) return CGAL_ZERO;
  return (CGAL_Sign) d.overlap_action();
}

inline CGAL_Comparison_result CGAL_compare(const CGAL_Interval_nt_advanced& d,
					   const CGAL_Interval_nt_advanced& e)
{
  if (e.sup < -d.inf) return CGAL_LARGER;
  if (-e.inf > d.sup) return CGAL_SMALLER;
  if ( (-d.inf==e.sup) && (-e.inf==d.sup) ) return CGAL_EQUAL;
  return (CGAL_Comparison_result) d.overlap_action();
}

inline CGAL_Interval_nt_advanced CGAL_abs (const CGAL_Interval_nt_advanced& d)
{
  if (d.inf <= 0) return d;
  if (d.sup <= 0) return -d;
  return CGAL_Interval_nt_advanced(0, CGAL_max(d.inf,d.sup));
}

inline CGAL_Interval_nt_advanced CGAL_min (const CGAL_Interval_nt_advanced& d,
					   const CGAL_Interval_nt_advanced& e)
{
  return CGAL_Interval_nt_advanced(-CGAL_max(d.inf, e.inf),
	  			    CGAL_min(d.sup, e.sup));
}

inline CGAL_Interval_nt_advanced CGAL_max (const CGAL_Interval_nt_advanced& d,
					   const CGAL_Interval_nt_advanced& e)
{
  return CGAL_Interval_nt_advanced(-CGAL_min(d.inf, e.inf),
	  			    CGAL_max(d.sup, e.sup));
}

ostream& operator<<(ostream& os, const CGAL_Interval_nt_advanced& d)
{ return os << "[" << d.lower_bound() << ";" << d.upper_bound() << "]"; }


// The non-advanced class.

struct CGAL_Interval_nt : public CGAL_Interval_nt_advanced
{
  friend CGAL_Interval_nt sqrt(const CGAL_Interval_nt&);
  friend CGAL_Interval_nt operator* (const double &, const CGAL_Interval_nt &);

public:

  // Constructors are identical.
  CGAL_Interval_nt()
      {}
  CGAL_Interval_nt(const double d)
      : CGAL_Interval_nt_advanced(d) {}
  CGAL_Interval_nt(const double a, const double b)
      : CGAL_Interval_nt_advanced(a,b) {}

  // This particular one needs to be redefined, a pitty...
  CGAL_Interval_nt operator-() const 
  { return CGAL_Interval_nt(-(sup), inf); }

  // The member functions that have to be protected against rounding mode.
  CGAL_Interval_nt operator+(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator-(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator*(const CGAL_Interval_nt& d) const ;
  CGAL_Interval_nt operator/(const CGAL_Interval_nt& d) const ;
  // For speed...
  CGAL_Interval_nt operator*(const double& d) const;
  // These have exactly the same code as the advanced class.
  // How can I avoid duplicating the code ?
  CGAL_Interval_nt& operator+=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator-=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator*=(const CGAL_Interval_nt& d) ;
  CGAL_Interval_nt& operator/=(const CGAL_Interval_nt& d) ;

private:
  // Private constructor for casts.
  CGAL_Interval_nt(const CGAL_Interval_nt_advanced &d)
      : CGAL_Interval_nt_advanced(d) {}
};

// Here we use the GNU extension of "Named return value".

#ifdef __GNUG__
#  define CGAL_NAMED_RETURN_VALUE_OPT_1 return tmp;
#  define CGAL_NAMED_RETURN_VALUE_OPT_2
#  define CGAL_NAMED_RETURN_VALUE_OPT_3
#else
#  define CGAL_NAMED_RETURN_VALUE_OPT_1
#  define CGAL_NAMED_RETURN_VALUE_OPT_2 CGAL_Interval_nt tmp;
#  define CGAL_NAMED_RETURN_VALUE_OPT_3 return tmp;
#endif

inline CGAL_Interval_nt CGAL_Interval_nt::operator+(const CGAL_Interval_nt& d)
    const CGAL_NAMED_RETURN_VALUE_OPT_1
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  tmp.inf = inf + d.inf;
  tmp.sup = sup + d.sup;
  CGAL_FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator-(const CGAL_Interval_nt& d)
    const CGAL_NAMED_RETURN_VALUE_OPT_1
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  tmp.inf = inf + d.sup;
  tmp.sup = sup + d.inf;
  CGAL_FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator*(const CGAL_Interval_nt& d)
    const
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator*(d) );
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt CGAL_Interval_nt::operator* (const double& d)
    const CGAL_NAMED_RETURN_VALUE_OPT_1
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  if (d>=0) {
      tmp.inf = inf*d;
      tmp.sup = sup*d;
  } else {
      tmp.inf = sup*(-d);
      tmp.sup = inf*(-d);
  }
  CGAL_FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline CGAL_Interval_nt operator*
  (const double& d, const CGAL_Interval_nt &t)
{ return t*d; }

inline CGAL_Interval_nt CGAL_Interval_nt::operator/(const CGAL_Interval_nt& d)
    const
{
  CGAL_FPU_set_rounding_to_infinity();
  CGAL_Interval_nt tmp ( CGAL_Interval_nt_advanced::operator/(d) );
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt &
  CGAL_Interval_nt::operator+=(const CGAL_Interval_nt& d)
{ return *this = *this + d; }

inline CGAL_Interval_nt &
  CGAL_Interval_nt::operator-=(const CGAL_Interval_nt& d)
{ return *this = *this - d; }

inline CGAL_Interval_nt &
  CGAL_Interval_nt::operator*=(const CGAL_Interval_nt& d)
{ return *this = *this * d; }

inline CGAL_Interval_nt &
  CGAL_Interval_nt::operator/=(const CGAL_Interval_nt& d)
{ return *this = *this / d; }

inline CGAL_Interval_nt sqrt(const CGAL_Interval_nt& d)
{
  CGAL_Interval_nt tmp = sqrt( (CGAL_Interval_nt_advanced) d);
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline CGAL_Interval_nt_advanced CGAL_to_interval_nt(const double &d)
{ return (CGAL_Interval_nt_advanced) d; }

// Finally we source the "cast" functions from other NTs, when necessary.

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

#endif // CGAL_INTERVAL_ARITHMETIC_H
