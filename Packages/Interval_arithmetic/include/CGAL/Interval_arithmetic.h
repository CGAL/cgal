// ============================================================================
//
// Copyright (c) 1998,1999 The CGAL Consortium
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


#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

// This file contains the description of the two classes:
// - Interval_nt_advanced  (do the FPU rounding mode changes yourself)
// - Interval_nt		("plug-in" version, derived from the other one)
//
// The differences are:
// - The second one is slower.
// - The first one supposes the rounding mode is set -> +infinity before
// nearly all operations, and might set it -> +infinity when leaving, whereas
// the second leaves the rounding -> nearest.

#include <iostream>
#include <CGAL/assertions.h>
#include <CGAL/IO/io_tags.h>		// For io_Operator().
#include <CGAL/number_type_tags.h>	// For number_type_tag()
#include <CGAL/double.h>	// For is_valid() and is_finite().
#include <CGAL/enum.h>  // Because we overload {sign,compare,abs,min,max}
#include <CGAL/number_utils.h> 		// For max and square<double>
#include <CGAL/Interval_arithmetic/_FPU.h>	// FPU rounding mode functions.
#include <CGAL/misc.h>			// For convert_to<>()

CGAL_BEGIN_NAMESPACE

struct Interval_nt_advanced
{
  typedef Interval_nt_advanced IA;
  struct unsafe_comparison{};			// Exception class.
  static const double min_double, max_double;	// Usefull constants.
  static const IA largest, smallest;

protected:
  double inf, sup;	// "inf" stores the OPPOSITE of the lower bound.
			// "sup" stores the upper bound of the interval.
private:
  int overlap_action() const
#ifndef CGAL_IA_NO_EXCEPTION
      throw (unsafe_comparison)
  { throw unsafe_comparison(); }
#else
  {
#if !defined(CGAL_IA_NO_WARNINGS) && !defined(CGAL_NO_WARNINGS)
     CGAL_warning_msg(false, " Comparison between overlapping intervals");
#endif
     return 0; // Return an arbitrary value.
  }
#endif // CGAL_IA_NO_EXCEPTION

public:
  friend IA	sqrt     (const IA &);
  friend IA	square	(const IA &);
  friend IA	abs (const IA &);
  friend IA	min (const IA &, const IA &);
  friend IA	max (const IA &, const IA &);
  friend IA	operator- (const double, const IA &);
  friend IA	operator/ (const double, const IA &);
  friend double to_double (const IA &);
  friend bool   is_valid  (const IA &);
  friend bool   is_finite (const IA &);
  friend Sign sign   (const IA &);
  friend Comparison_result compare (const IA &, const IA &);

  // The constructors.
  Interval_nt_advanced()
#ifndef CGAL_NO_ASSERTIONS
      : inf(-1), sup(-1) // Buggy interval to detect use before definition.
#endif
      {}

  Interval_nt_advanced(const double d)
      : inf(-d), sup(d) {}

  Interval_nt_advanced(const double i, const double s)
      : inf(-i), sup(s)
      { CGAL_assertion_msg(i<=s," Variable used before being initialized ?"); }

#if 1 // Need to be redefined for non advanced ?
  // The copy constructors/assignment: useless.
  // The default ones are ok, but these are faster...
  Interval_nt_advanced(const IA & d)
      : inf(d.inf), sup(d.sup) {}

  IA & operator=(const IA & d)
  { inf = d.inf; sup = d.sup; return *this; }
#endif

  // The operators.
  IA  operator+(const IA & d) const
  {
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
      return IA (-(inf + d.inf), sup + d.sup);
  }
  // { return IA  (d) += *this; }

  IA  operator-(const IA & d) const
  {
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
      return IA (-(inf + d.sup), sup + d.inf);
  }

  // Those 2 ones could be made not inlined.
  IA  operator*(const IA & d) const;
  IA  operator/(const IA & d) const;

  IA  operator-() const { return IA (-(sup), inf); }

  IA & operator+=(const IA & d);
  IA & operator-=(const IA & d);
  IA & operator*=(const IA & d);
  IA & operator/=(const IA & d);

  // For speed...
  IA  operator+(const double d) const {  return IA(-(inf-d), sup+d); };
  IA  operator-(const double d) const {  return IA(-(inf+d), sup-d); };
  IA  operator*(const double d) const;
  IA  operator/(const double d) const;

  bool operator<(const IA & d) const
  {
    if (sup  < -d.inf) return true;
    if (-inf >= d.sup) return false;
    return overlap_action();
  }

  bool operator<=(const IA & d) const
  {
    if (sup <= -d.inf) return true;
    if (-inf >  d.sup) return false;
    return overlap_action();
  }
  
  bool operator==(const IA & d) const
  {
    if ((-d.inf >  sup) || (d.sup  < -inf)) return false;
    if ((-d.inf == sup) && (d.sup == -inf)) return true;
    return overlap_action();
  }

  bool operator> (const IA & d) const { return  (d <  *this); }
  bool operator>=(const IA & d) const { return  (d <= *this); }
  bool operator!=(const IA & d) const { return !(d == *this); }

  bool is_same(const IA & d) const
  { return (inf == d.inf) && (sup == d.sup); }

  bool overlap(const IA & d) const
  { return !((-d.inf > sup) || (d.sup < -inf)); }

  double lower_bound() const { return -inf; }
  double upper_bound() const { return sup; }
};

// Usefull constants.

// MIN_DOUBLE (subnormal)
const double Interval_nt_advanced::min_double =5e-324;

// MAX_DOUBLE
const double Interval_nt_advanced::max_double =1.7976931348623157081e+308;

// Smallest interval strictly around zero.
const Interval_nt_advanced Interval_nt_advanced::smallest
(-Interval_nt_advanced::min_double,
  Interval_nt_advanced::min_double);

// [-inf;+inf]
const Interval_nt_advanced Interval_nt_advanced::largest
(-HUGE_VAL, HUGE_VAL);


inline
Interval_nt_advanced
Interval_nt_advanced::operator* (const Interval_nt_advanced & d) const
{
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
  if (inf<=0)					// this>=0
  {
      /* d>=0     return IA (-((-inf)*d.inf),    sup*d.sup);
       * d<=0     return IA (-(   sup*d.inf), (-inf)*d.sup);
       * 0 \in d  return IA (-(   sup*d.inf),    sup*d.sup);
       */
    double a = -inf, b = sup;
    if (d.inf > 0)
    {
	a=b;
	if (d.sup < 0)
	    b=-inf;
    }
    return IA (-(a*d.inf), b*d.sup);
  }
  else if (sup<=0)				// this<=0
  {
      /* d>=0     return IA (-(   inf*d.sup), (-sup)*d.inf);
       * d<=0     return IA (-((-sup)*d.sup),    inf*d.inf);
       * 0 \in d  return IA (-(   inf*d.sup),    inf*d.inf);
       */
    double a = -sup, b = inf;
    if (d.inf > 0)
    {
	a=b;
	if (d.sup < 0)
	    b=-sup;
    }
    return IA (-(b*d.sup), a*d.inf);
  }
  else						// 0 \in [inf;sup]
  {
    if (d.inf<=0)				// d>=0
      return IA (-(inf*d.sup), sup*d.sup);
    if (d.sup<=0)				// d<=0
      return IA (-(sup*d.inf), inf*d.inf);
        					// 0 \in d
    double tmp1 = inf*d.sup;
    double tmp2 = sup*d.inf;
    double tmp3 = inf*d.inf;
    double tmp4 = sup*d.sup;
    return IA (-max(tmp1,tmp2), max(tmp3,tmp4));
  };
}

inline
Interval_nt_advanced
Interval_nt_advanced::operator* (const double d) const
{
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
#if 0
      double a = inf*fabs(d);
      double b = sup*fabs(d);
      if (d>=0) return IA(-a,b);
      return IA(-b,a);
#else
  if (d>=0) return IA (-(inf*d), sup*d);
  return IA (-(sup*(-d)), inf*(-d));
#endif
}

inline
Interval_nt_advanced
Interval_nt_advanced::operator/ (const double d) const
{
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
  if (d>0) return IA (-(inf/d), sup/d);
  if (d<0) return IA (-(sup/(-d)), inf/(-d));
  return largest;
}

inline
bool
operator< (const double d, const Interval_nt_advanced & t)
{ return t>d; }

inline
bool
operator<= (const double d, const Interval_nt_advanced & t)
{ return t>=d; }

inline
bool
operator> (const double d, const Interval_nt_advanced & t)
{ return t<d; }

inline
bool
operator>= (const double d, const Interval_nt_advanced & t)
{ return t<=d; }

inline
bool
operator== (const double d, const Interval_nt_advanced & t)
{ return t==d; }

inline
bool
operator!= (const double d, const Interval_nt_advanced & t)
{ return t!=d; }

inline
Interval_nt_advanced
operator+ (const double d, const Interval_nt_advanced & t)
{ return t+d; }

inline
Interval_nt_advanced
operator- (const double d, const Interval_nt_advanced & t)
{ return Interval_nt_advanced(-(t.sup-d), t.inf+d); }
// { return -(t-d); }

inline
Interval_nt_advanced
operator* (const double d, const Interval_nt_advanced & t)
{ return t*d; }

inline
Interval_nt_advanced
operator/ (const double t, const Interval_nt_advanced & d)
{
  typedef Interval_nt_advanced IA;
  if (d.inf<0)				// d>0
  {
    if (t>=0) return IA(-((-t)/d.sup), (-t)/d.inf);
              return IA(-(t/d.inf),    t/d.sup);
  }
  if (d.sup<0)				// d<0
  {
    if (t>=0) return IA(-(t/d.inf),    t/d.sup);
              return IA(-((-t)/d.sup), (-t)/d.inf);
  }
  return IA::largest;			// 0 \in d
}

inline
Interval_nt_advanced
Interval_nt_advanced::operator/ (const Interval_nt_advanced & d) const
{
#ifdef CGAL_IA_DEBUG
      CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
  if (d.inf<0)				// d>0
  {
      /* this>=0	return IA (-(inf/d.sup), sup/(-d.inf));
       * this<=0	return IA (-(inf/(-d.inf)), sup/d.sup);
       * 0 \in this	return IA (-(inf/(-d.inf)), sup/(-d.inf));
       */
    double a = d.sup, b = -d.inf;
    if (inf>0)
    {
	a=b;
	if (sup<0)
	    b=d.sup;
    };
    return IA(-(inf/a), sup/b);
  }
  else if (d.sup<0)			// d<0
  {
      /* this>=0	return IA (-(sup/(-d.sup)), inf/d.inf);
       * this<=0	return IA (-(sup/d.inf),    inf/(-d.sup));
       * 0 \in this	return IA (-(sup/(-d.sup)), inf/(-d.sup));
       */
    double a = -d.sup, b = d.inf;
    if (inf>0)
    {
	a=b;
	if (sup<0)
	    b=-d.sup;
    };
    return IA(-(sup/a), inf/b);
  }
  else					// 0 \in d
    return largest; // IA (-HUGE_VAL, HUGE_VAL);
	   // We could do slightly better -> [0;HUGE_VAL] when d.sup==0,
	   // but is this worth ?
}

inline
Interval_nt_advanced &
Interval_nt_advanced::operator+= (const Interval_nt_advanced & d)
{ return *this = *this + d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator-= (const Interval_nt_advanced & d)
{ return *this = *this - d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator*= (const Interval_nt_advanced & d)
{ return *this = *this * d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator/= (const Interval_nt_advanced & d)
{ return *this = *this / d; }

inline
Interval_nt_advanced
sqrt (const Interval_nt_advanced & d)
{
  FPU_set_rounding_to_minus_infinity();
  Interval_nt_advanced tmp;
  // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
  // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
  // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
  tmp.inf = (d.inf<0) ? -sqrt(-d.inf) : 0;
  FPU_set_rounding_to_infinity();
  tmp.sup = sqrt(d.sup);
  return tmp;
}

inline
Interval_nt_advanced
square (const Interval_nt_advanced & d)
{
    // The first one is slightly slower, but produces shorter code.
#if 0
  double a = d.inf*fabs(d.inf);
  double b = d.sup*fabs(d.sup);
  if (d.inf <= 0) return Interval_nt_advanced(-a,b);
  if (d.sup <= 0) return Interval_nt_advanced(-b,a);
  return Interval_nt_advanced(0, max(a,b));
#else
  if (d.inf<=0) return Interval_nt_advanced(-(d.inf*-d.inf), d.sup*d.sup);
  if (d.sup<=0) return Interval_nt_advanced(-(d.sup*-d.sup), d.inf*d.inf);
  return Interval_nt_advanced(0.0, square(max(d.inf, d.sup)));
#endif
}

inline
double
to_double (const Interval_nt_advanced & d)
{ return (d.sup-d.inf)*.5; }

inline
bool
is_valid (const Interval_nt_advanced & d)
{ return is_valid(d.inf) && is_valid(d.sup) && (-d.inf <= d.sup); }

inline
bool
is_finite (const Interval_nt_advanced & d)
{ return is_finite(d.inf) && is_finite(d.sup); }

inline
Sign
sign (const Interval_nt_advanced & d)
{
  if (d.inf < 0) return POSITIVE;
  if (d.sup < 0) return NEGATIVE;
  if (-d.inf == d.sup) return ZERO;
  return Sign (d.overlap_action());
}

inline
Comparison_result
compare (const Interval_nt_advanced & d,
	const Interval_nt_advanced & e)
{
  if (-d.inf > e.sup) return LARGER;
  if (-e.inf > d.sup) return SMALLER;
  if ( (-e.inf == d.sup) && (-d.inf == e.sup) ) return EQUAL;
  return Comparison_result (d.overlap_action());
}

inline
Interval_nt_advanced
abs (const Interval_nt_advanced & d)
{
  if (d.inf <= 0) return d;
  if (d.sup <= 0) return -d;
  return Interval_nt_advanced(0, max(d.inf,d.sup));
}

inline
Interval_nt_advanced
min (const Interval_nt_advanced & d,
	const Interval_nt_advanced & e)
{
  return Interval_nt_advanced(-max(d.inf, e.inf),
	  			    min(d.sup, e.sup));
}

inline
Interval_nt_advanced
max (const Interval_nt_advanced & d,
	const Interval_nt_advanced & e)
{
  return Interval_nt_advanced(-min(d.inf, e.inf),
	  			    max(d.sup, e.sup));
}

inline
ostream &
operator<< (ostream & os, const Interval_nt_advanced & d)
{ return os << "[" << d.lower_bound() << ";" << d.upper_bound() << "]"; }

inline
istream &
operator>> (istream & is, Interval_nt_advanced & ia)
{
    double d;
    is >> d;
    ia = d;
    return is;
}


// The non-advanced class.

struct Interval_nt : public Interval_nt_advanced
{
  typedef Interval_nt IA;

  // Constructors are identical.
  Interval_nt()
      {}
  Interval_nt(const double d)
      : Interval_nt_advanced(d) {}
  Interval_nt(const double a, const double b)
      : Interval_nt_advanced(a,b) {}

  // Private constructor for casts. (remade public)
  Interval_nt(const Interval_nt_advanced &d)
      : Interval_nt_advanced(d) {}

  friend IA	sqrt	(const IA &);
  friend IA	square	(const IA &);
  // How to share those definitions ?
  friend IA     abs (const IA & d)
    { return abs((Interval_nt_advanced) d); }
  friend IA     min (const IA & d, const IA & e)
    { return min((Interval_nt_advanced) d,
	              (Interval_nt_advanced) e);
    }
  friend IA     max (const IA & d, const IA & e)
    { return max((Interval_nt_advanced) d,
	              (Interval_nt_advanced) e);
    }
  // friend IA     operator* (const double &, const IA &);
  friend double to_double (const IA & d)
    { return to_double((Interval_nt_advanced) d); }
  friend bool   is_valid  (const IA & d)
    { return is_valid((Interval_nt_advanced) d); }
  friend bool   is_finite (const IA & d)
    { return is_finite((Interval_nt_advanced) d); }
  friend Sign sign   (const IA & d)
    { return sign((Interval_nt_advanced) d); }
  friend Comparison_result compare (const IA & d, const IA & e)
    { return compare((Interval_nt_advanced) d,
	                  (Interval_nt_advanced) e);
    }

  // This particular one needs to be redefined, a pitty...
  IA operator-() const 
  { return IA(-(sup), inf); }

  // The member functions that have to be protected against rounding mode.
  IA operator+(const IA & d) const ;
  IA operator-(const IA & d) const ;
  IA operator*(const IA & d) const ;
  IA operator/(const IA & d) const ;
  // For speed...
  IA operator*(const double d) const;
  // These have exactly the same code as the advanced class.
  // How can I avoid duplicating the code ?
  IA & operator+=(const IA & d) ;
  IA & operator-=(const IA & d) ;
  IA & operator*=(const IA & d) ;
  IA & operator/=(const IA & d) ;
};

// Here we use the GNU extension of "Named return value".

#ifdef __GNUG__
#  define CGAL_NAMED_RETURN_VALUE_OPT_1 return tmp;
#  define CGAL_NAMED_RETURN_VALUE_OPT_2
#  define CGAL_NAMED_RETURN_VALUE_OPT_3
#else
#  define CGAL_NAMED_RETURN_VALUE_OPT_1
#  define CGAL_NAMED_RETURN_VALUE_OPT_2 Interval_nt tmp;
#  define CGAL_NAMED_RETURN_VALUE_OPT_3 return tmp;
#endif

inline
Interval_nt
Interval_nt::operator+ (const Interval_nt & d) const
CGAL_NAMED_RETURN_VALUE_OPT_1
{
  FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  tmp.inf = inf + d.inf;
  tmp.sup = sup + d.sup;
  FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline
Interval_nt
Interval_nt::operator- (const Interval_nt & d) const
CGAL_NAMED_RETURN_VALUE_OPT_1
{
  FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  tmp.inf = inf + d.sup;
  tmp.sup = sup + d.inf;
  FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline
Interval_nt
Interval_nt::operator* (const Interval_nt & d) const
{
  FPU_set_rounding_to_infinity();
  Interval_nt tmp ( Interval_nt_advanced::operator*(d) );
  FPU_set_rounding_to_nearest();
  return tmp;
}

inline
Interval_nt
Interval_nt::operator* (const double d) const
CGAL_NAMED_RETURN_VALUE_OPT_1
{
  FPU_set_rounding_to_infinity();
  CGAL_NAMED_RETURN_VALUE_OPT_2
  if (d>=0) {
      tmp.inf = inf*d;
      tmp.sup = sup*d;
  } else {
      tmp.inf = sup*(-d);
      tmp.sup = inf*(-d);
  }
  FPU_set_rounding_to_nearest();
  CGAL_NAMED_RETURN_VALUE_OPT_3
}

inline
Interval_nt
operator* (const double d, const Interval_nt & t)
{ return t*d; }

inline
Interval_nt
Interval_nt::operator/ (const Interval_nt & d) const
{
  FPU_set_rounding_to_infinity();
  Interval_nt tmp ( Interval_nt_advanced::operator/(d) );
  FPU_set_rounding_to_nearest();
  return tmp;
}

inline
Interval_nt &
Interval_nt::operator+= (const Interval_nt & d)
{ return *this = *this + d; }

inline
Interval_nt &
Interval_nt::operator-= (const Interval_nt & d)
{ return *this = *this - d; }

inline
Interval_nt &
Interval_nt::operator*= (const Interval_nt & d)
{ return *this = *this * d; }

inline
Interval_nt &
Interval_nt::operator/= (const Interval_nt & d)
{ return *this = *this / d; }

inline
Interval_nt
sqrt (const Interval_nt & d)
{
  Interval_nt tmp = sqrt( (Interval_nt_advanced) d);
  FPU_set_rounding_to_nearest();
  return tmp;
}

inline
Interval_nt
square (const Interval_nt & d)
{
  FPU_set_rounding_to_infinity();
  Interval_nt tmp = square( (Interval_nt_advanced) d);
  FPU_set_rounding_to_nearest();
  return tmp;
}


// The undocumented Tag things...

inline
io_Operator
io_tag (const Interval_nt_advanced &)
{ return io_Operator(); }

inline
Number_tag
number_type_tag (Interval_nt_advanced)
{ return Number_tag(); }



// Finally we deal with the convert_to<Interval_nt_advanced>(NT)
// functions from other NTs, when necessary.
// convert_to<Interval_nt>() is templated below.
//
// For the builtin types (well, all those that can be casted to double
// exactly), the template in misc.h is enough.

CGAL_END_NAMESPACE

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

#ifdef CGAL_FIXED_PRECISION_NT_H
#include <CGAL/Interval_arithmetic/IA_Fixed.h>
#endif

#ifdef CGAL_QUOTIENT_H
#include <CGAL/Interval_arithmetic/IA_Quotient.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class FT>
inline
Interval_nt
convert_to (const FT & z)
{
    FPU_set_rounding_to_infinity();
    Interval_nt tmp(convert_to<Interval_nt_advanced>(z));
    FPU_set_rounding_to_nearest();
    return tmp;
}

CGAL_END_NAMESPACE

#endif // CGAL_INTERVAL_ARITHMETIC_H
