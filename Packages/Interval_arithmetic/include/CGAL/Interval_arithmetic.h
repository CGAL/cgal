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
//
// The differences are:
// - The second one is slower.
// - The first one supposes the rounding mode is set -> +infinity before
// nearly all operations, and might set it -> +infinity when leaving, whereas
// the second leaves the rounding -> nearest.

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

#include <iostream.h>
#include <CGAL/assertions.h>
#include <CGAL/double.h>	// For CGAL_is_valid() and CGAL_is_finite().
#include <CGAL/enum.h>  // Because we overload CGAL_{sign,compare,abs,min,max}
#include <CGAL/number_utils.h>  // For CGAL_max(double, double).
#include <CGAL/Interval_arithmetic/_FPU.h>	// FPU rounding mode functions.

struct CGAL_Interval_nt_advanced
{
  typedef CGAL_Interval_nt_advanced IA;
  struct unsafe_comparison{};  // Exception class.
  // Usefull constants.
  const static double min_double = 5e-324; // MIN_DOUBLE (subnormal)
  const static double max_double = 1.7976931348623157081e+308; // MAX_DOUBLE

protected:
  double inf, sup;	// "inf" stores the OPPOSITE of the lower bound.
			// "sup" stores the upper bound of the interval.
private:
  int overlap_action() const throw (unsafe_comparison)
  { throw unsafe_comparison(); }

public:
  friend IA	sqrt     (const IA &);
  friend IA	CGAL_abs (const IA &);
  friend IA	CGAL_min (const IA &, const IA &);
  friend IA	CGAL_max (const IA &, const IA &);
  // friend IA	operator* (const double, const IA &);
  friend double CGAL_to_double (const IA &);
  friend bool   CGAL_is_valid  (const IA &);
  friend bool   CGAL_is_finite (const IA &);
  friend CGAL_Sign CGAL_sign   (const IA &);
  friend CGAL_Comparison_result CGAL_compare (const IA &, const IA &);

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
  CGAL_Interval_nt_advanced(const IA & d)
      : inf(d.inf), sup(d.sup) {}

  IA & operator=(const IA & d)
  { inf = d.inf; sup = d.sup; return *this; }
#endif

  // The operators.
  IA  operator+(const IA & d) const
  { return IA (-(inf + d.inf), sup + d.sup); }
  // { return IA  (d) += *this; }

  IA  operator-(const IA & d) const
  { return IA (-(inf + d.sup), sup + d.inf); }

  // Those 2 ones could be made not inlined.
  IA  operator*(const IA & d) const;
  IA  operator/(const IA & d) const;

  IA  operator-() const
  { return IA (-(sup), inf); }

  IA & operator+=(const IA & d);
  IA & operator-=(const IA & d);
  IA & operator*=(const IA & d);
  IA & operator/=(const IA & d);

  // For speed...
  IA  operator*(const double d) const;

  bool operator<(const IA & d) const
  {
    if (sup  < -d.inf) return true;
    if (-inf >= d.sup) return false;
    return overlap_action();
  }

  bool operator>(const IA & d) const
  { return (d < *this); }

  bool operator<=(const IA & d) const
  {
    if (sup <= -d.inf) return true;
    if (-inf >  d.sup) return false;
    return overlap_action();
  }

  bool operator>=(const IA & d) const
  { return (d <= *this); }

  bool operator==(const IA & d) const
  {
    if ((-d.inf >  sup) || (d.sup  < -inf)) return false;
    if ((-d.inf == sup) && (d.sup == -inf)) return true;
    return overlap_action();
  }

  bool operator!=(const IA & d) const
  { return !(d == *this); }
  
  bool is_same(const IA & d) const
  { return (inf == d.inf) && (sup == d.sup); }

  double lower_bound() const { return -inf; }
  double upper_bound() const { return sup; }
};

// Two usefull constants.

static const CGAL_Interval_nt_advanced CGAL_Interval_smallest
(-CGAL_Interval_nt_advanced::min_double, CGAL_Interval_nt_advanced::min_double);

static const CGAL_Interval_nt_advanced CGAL_Interval_largest
(-HUGE_VAL, HUGE_VAL);


inline CGAL_Interval_nt_advanced CGAL_Interval_nt_advanced::operator*
  (const CGAL_Interval_nt_advanced& d) const
{
  if (inf<=0)					/* this>=0 */
  {
    if (d.inf<=0)				/* d>=0 */
      return IA (-((-inf)*d.inf), sup*d.sup);
    else if (d.sup<=0)				/* d<=0 */
      return IA (-(sup*d.inf), (-inf)*d.sup);
    else					/* 0 \in d */
      return IA (-(sup*d.inf), sup*d.sup);
  }
  else if (sup<=0)				/* this<=0 */
  {
    if (d.inf<=0)				/* d>=0 */
      return IA (-(inf*d.sup), sup*(-d.inf));
    else if (d.sup<=0)				/* d<=0 */
      return IA (-((-sup)*d.sup), inf*d.inf);
    else					/* 0 \in d */
      return IA (-(inf*d.sup), inf*d.inf);
  }
  else						/* 0 \in [inf;sup] */
  {
    if (d.inf<=0)				/* d>=0 */
      return IA (-(inf*d.sup), sup*d.sup);
    else if (d.sup<=0)				/* d<=0 */
      return IA (-(sup*d.inf), inf*d.inf);
    else					/* 0 \in d */
    {
      double tmp1, tmp2, tmp3, tmp4;
      tmp1 = inf*d.sup;
      tmp2 = sup*d.inf;
      tmp3 = inf*d.inf;
      tmp4 = sup*d.sup;
      return IA (-CGAL_max(tmp1,tmp2), CGAL_max(tmp3,tmp4));
    };
  };
}

inline CGAL_Interval_nt_advanced CGAL_Interval_nt_advanced::operator*
  (const double d) const
{
  if (d>=0)	return IA (-(inf*d), sup*d);
  else		return IA (-(sup*(-d)), inf*(-d));
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
      return IA (-(inf/d.sup), sup/(-d.inf));
    else if (sup<=0.0)				/* this<=0 */
      return IA (-(inf/(-d.inf)), sup/d.sup);
    else					/* 0 \in this */
      return IA (-(inf/(-d.inf)), sup/(-d.inf));
  }
  else if (d.sup<0.0)				/* d<0 */
  {
    if (inf<=0.0)				/* this>=0 */
      return IA (-(sup/(-d.sup)), inf/d.inf);
    else if (sup<=0.0)				/* this<=0 */
      return IA (-(sup/d.inf), inf/(-d.sup));
    else					/* 0 \in this */
      return IA (-(sup/(-d.sup)), inf/(-d.sup));
  }
  else						/* 0 \in [d.inf;d.sup] */
    return CGAL_Interval_largest; // IA (-HUGE_VAL, HUGE_VAL);
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
  // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
  // sqrt([-a,+b]) => [       0;sqrt(+b)] => we assume roundoff error.
  // sqrt([-a,-b]) => [sqrt(-a);sqrt(-b)] => we assume user bug.
  if ((d.inf<0) || (d.sup<0))
    tmp.inf = - sqrt(-d.inf);
  else
    tmp.inf = 0;
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
private:

  // For this one to be overriden, it must be virtual...
#if 0
  int overlap_action() const
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
     CGAL_warning_msg(false, " Comparison between overlapping intervals");
#endif
     return 0; // Return a "random" value.
  }
#endif

public:

  // Constructors are identical.
  CGAL_Interval_nt()
      {}
  CGAL_Interval_nt(const double d)
      : CGAL_Interval_nt_advanced(d) {}
  CGAL_Interval_nt(const double a, const double b)
      : CGAL_Interval_nt_advanced(a,b) {}

  // Private constructor for casts. (remade public)
  CGAL_Interval_nt(const CGAL_Interval_nt_advanced &d)
      : CGAL_Interval_nt_advanced(d) {}

  typedef CGAL_Interval_nt IA;

  friend IA     sqrt     (const IA &);
  // How to share those definitions ?
  friend IA     CGAL_abs (const IA & d)
    { return CGAL_abs((CGAL_Interval_nt_advanced) d); }
  friend IA     CGAL_min (const IA & d, const IA & e)
    { return CGAL_min((CGAL_Interval_nt_advanced) d,
	              (CGAL_Interval_nt_advanced) e);
    }
  friend IA     CGAL_max (const IA & d, const IA & e)
    { return CGAL_max((CGAL_Interval_nt_advanced) d,
	              (CGAL_Interval_nt_advanced) e);
    }
  // friend IA     operator* (const double &, const IA &);
  friend double CGAL_to_double (const IA &);
  friend bool   CGAL_is_valid  (const IA & d)
    { return CGAL_is_valid((CGAL_Interval_nt_advanced) d); }
  friend bool   CGAL_is_finite (const IA & d)
    { return CGAL_is_finite((CGAL_Interval_nt_advanced) d); }
  friend CGAL_Sign CGAL_sign   (const IA & d)
    { return CGAL_sign((CGAL_Interval_nt_advanced) d); }
  friend CGAL_Comparison_result CGAL_compare (const IA &d, const IA &e)
    { return CGAL_compare((CGAL_Interval_nt_advanced) d,
	                  (CGAL_Interval_nt_advanced) e);
    }

  // This particular one needs to be redefined, a pitty...
  IA operator-() const 
  { return IA(-(sup), inf); }

  // The member functions that have to be protected against rounding mode.
  IA operator+(const IA& d) const ;
  IA operator-(const IA& d) const ;
  IA operator*(const IA& d) const ;
  IA operator/(const IA& d) const ;
  // For speed...
  IA operator*(const double d) const;
  // These have exactly the same code as the advanced class.
  // How can I avoid duplicating the code ?
  IA& operator+=(const IA& d) ;
  IA& operator-=(const IA& d) ;
  IA& operator*=(const IA& d) ;
  IA& operator/=(const IA& d) ;
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

inline CGAL_Interval_nt CGAL_Interval_nt::operator* (const double d)
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
  (const double d, const CGAL_Interval_nt &t)
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

// Finally we source the CGAL_to_Interval_nt_advanced() "cast" functions from
// other NTs, when necessary.
// CGAL_to_Interval_nt() is templated below, only specialized for double.

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

inline CGAL_Interval_nt_advanced CGAL_to_Interval_nt_advanced(const double d)
{ return CGAL_Interval_nt_advanced(d); }

template <class FT>
inline CGAL_Interval_nt CGAL_to_Interval_nt(const FT &z)
{
    CGAL_FPU_set_rounding_to_infinity();
    CGAL_Interval_nt tmp(CGAL_to_Interval_nt_advanced(z));
    CGAL_FPU_set_rounding_to_nearest();
    return tmp;
}

inline CGAL_Interval_nt CGAL_to_Interval_nt(const double d)
{ return CGAL_Interval_nt(d); }

#endif // CGAL_INTERVAL_ARITHMETIC_H
