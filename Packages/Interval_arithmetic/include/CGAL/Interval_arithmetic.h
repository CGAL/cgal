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
// revision      : 1.5
// revision_date : 26 February 1998
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Herve.Bronnimann@sophia.inria.fr>)
//
// ============================================================================

// This file contains the description of the two classes:
// - CGAL_Interval_nt_advanced  (do the FPU tricks yourself)
// - CGAL_Interval_nt		(derived from the other one)

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

#include <iostream.h>
#include <assert.h>
#include <CGAL/double.h>        // For CGAL_is_valid() and CGAL_is_finite().
#include <CGAL/_FPU.h>          // FPU rounding mode functions.


class CGAL_Interval_nt_advanced
{
  friend CGAL_Interval_nt_advanced sqrt(const CGAL_Interval_nt_advanced&);
  friend double CGAL_to_double(const CGAL_Interval_nt_advanced& d);

public:

  inline CGAL_Interval_nt_advanced(double i, double s)
  {
#ifndef CGAL_NO_PRECONDITIONS
    assert(i<=s);
#endif
    inf = -i; sup = s;
  }

  inline CGAL_Interval_nt_advanced(double d = 0.0)
  {
    inf = -d; sup = d;
  }

  inline CGAL_Interval_nt_advanced&
         operator=(const CGAL_Interval_nt_advanced& d)
  {
    inf = d.inf; sup = d.sup;
    return *this;
  }

  inline CGAL_Interval_nt_advanced
         operator+(const CGAL_Interval_nt_advanced& d) const
  {
    CGAL_Interval_nt_advanced tmp;
    tmp.inf = inf + d.inf;
    tmp.sup = sup + d.sup;
    return tmp;
  }

  inline CGAL_Interval_nt_advanced
         operator-(const CGAL_Interval_nt_advanced& d) const
  {
    CGAL_Interval_nt_advanced tmp;
    tmp.inf = inf + d.sup;
    tmp.sup = sup + d.inf;
    return tmp;
  }

  inline CGAL_Interval_nt_advanced
         operator*(const CGAL_Interval_nt_advanced& d) const
  {
    CGAL_Interval_nt_advanced tmp;
    if (inf<=0)                                   /* this>=0 */
    {
      if (d.inf<=0)                                 /* d>=0 */
      {
        tmp.inf = (-inf)*d.inf;
        tmp.sup = sup*d.sup;
      }
      else if (d.sup<=0)                            /* d<=0 */
      {
        tmp.inf = sup*d.inf;
        tmp.sup = (-inf)*d.sup;
      }
      else                                        /* 0 \in d */
      {
        tmp.inf = sup*d.inf;
        tmp.sup = sup*d.sup;
      };
    }
    else if (sup<=0)                              /* this<=0 */
    {
      if (d.inf<=0)                                 /* d>=0 */
      {
        tmp.inf = inf*d.sup;
        tmp.sup = sup*(-d.inf);
      }
      else if (d.sup<=0)                            /* d<=0 */
      {
        tmp.inf = (-sup)*d.sup;
        tmp.sup = inf*d.inf;
      }
      else                                        /* 0 \in d */
      {
        tmp.inf = inf*d.sup;
        tmp.sup = inf*d.inf;
      };
    }
    else                                          /* 0 \in [inf;sup] */
    {
      if (d.inf<=0)                                 /* d>=0 */
      {
        tmp.inf = inf*d.sup;
        tmp.sup = sup*d.sup;
      }
      else if (d.sup<=0)                            /* d<=0 */
      {
        tmp.inf = sup*d.inf;
        tmp.sup = inf*d.inf;
      }
      else                                        /* 0 \in d */
      {
        double tmp1, tmp2;
        tmp1 = inf*d.sup;
        tmp2 = sup*d.inf;
        tmp.inf = (tmp1>tmp2) ? tmp1 : tmp2;
        tmp1 = inf*d.inf;
        tmp2 = sup*d.sup;
        tmp.sup = (tmp1>tmp2) ? tmp1 : tmp2;
      };
    };
    return tmp;
  }

  inline CGAL_Interval_nt_advanced
         operator/(const CGAL_Interval_nt_advanced& d) const
  {
    CGAL_Interval_nt_advanced tmp;
    if (d.inf<0.0)                            /* d>0 */
    {
      if (inf<=0.0)                         /* this>=0 */
      {
        tmp.inf = inf/d.sup;
        tmp.sup = sup/(-d.inf);
      }
      else if (sup<=0.0)                    /* this<=0 */
      {
        tmp.inf = inf/(-d.inf);
        tmp.sup = sup/d.sup;
      }
      else                                /* 0 \in this */
      {
        tmp.inf = inf/(-d.inf);
        tmp.sup = sup/(-d.inf);
      };
    }
    else if (d.sup<0.0)                       /* d<0 */
    {
      if (inf<=0.0)                         /* this>=0 */
      {
        tmp.inf = sup/(-d.sup);
        tmp.sup = inf/d.inf;
      }
      else if (sup<=0.0)                    /* b<=0 */
      {
        tmp.inf = sup/d.inf;
        tmp.sup = inf/(-d.sup);
      }
      else                                /* 0 \in this */
      {
        tmp.inf = sup/(-d.sup);
        tmp.sup = inf/(-d.sup);
      };
    }
    else                                  /* 0 \in [d.inf;d.sup] */
      tmp.inf = tmp.sup = HUGE_VAL;           /* ]-oo;+oo[. */
    return tmp;
  }

  inline CGAL_Interval_nt_advanced operator-() const
  {
    CGAL_Interval_nt_advanced tmp;
    tmp.inf = sup;
    tmp.sup = inf;
    return tmp;
  }

  inline bool is_valid() const
  {
    return CGAL_is_valid(inf) && CGAL_is_valid(sup);
  }

  inline bool is_finite() const
  {
    return CGAL_is_finite(inf) && CGAL_is_finite(sup);
  }

  inline bool operator==(const CGAL_Interval_nt_advanced& d) const
  {
    return !(*this != d);
  }

  inline bool operator!=(const CGAL_Interval_nt_advanced& d) const
  {
    return (*this < d) || (d < *this);
  }

  inline bool operator<(const CGAL_Interval_nt_advanced& d) const
  {
    return (sup < -d.inf);
  }

  inline bool operator>(const CGAL_Interval_nt_advanced& d) const
  {
    return (d.sup < -inf);
  }

  inline bool operator<=(const CGAL_Interval_nt_advanced& d) const
  {
    return !(*this > d);
  }

  inline bool operator>=(const CGAL_Interval_nt_advanced& d) const
  {
    return !(*this < d);
  }

  inline double lower_bound() const
  {
    return -inf;
  }

  inline double upper_bound() const
  {
    return sup;
  }

protected:
        // "inf" stores the __opposite__ of the inferior bound.
        // "sup" stores the upper bound of the interval.
  double inf, sup;
};


class CGAL_Interval_nt : public CGAL_Interval_nt_advanced
{
  friend CGAL_Interval_nt sqrt(const CGAL_Interval_nt&);

public:

  inline CGAL_Interval_nt(double d = 0.0) : CGAL_Interval_nt_advanced(d) {}
  inline CGAL_Interval_nt(double a,double b) : CGAL_Interval_nt_advanced(a,b) {}

  inline CGAL_Interval_nt(const CGAL_Interval_nt_advanced &d)
    : CGAL_Interval_nt_advanced(d) {}

  inline CGAL_Interval_nt operator+(const CGAL_Interval_nt& d) const
  {
    CGAL_Interval_nt tmp;
    CGAL_FPU_set_rounding_to_infinity();
    tmp = ((CGAL_Interval_nt_advanced) *this) + (CGAL_Interval_nt_advanced) d;
    CGAL_FPU_set_rounding_to_nearest();
    return tmp;
  }

  inline CGAL_Interval_nt operator-(const CGAL_Interval_nt& d) const
  {
    CGAL_Interval_nt tmp;
    CGAL_FPU_set_rounding_to_infinity();
    tmp = ((CGAL_Interval_nt_advanced) *this) - (CGAL_Interval_nt_advanced) d;
    CGAL_FPU_set_rounding_to_nearest();
    return tmp;
  }

  inline CGAL_Interval_nt operator*(const CGAL_Interval_nt& d) const
  {
    CGAL_Interval_nt tmp;
    CGAL_FPU_set_rounding_to_infinity();
    tmp = ((CGAL_Interval_nt_advanced) *this) * (CGAL_Interval_nt_advanced) d;
    CGAL_FPU_set_rounding_to_nearest();
    return tmp;
  }

  inline CGAL_Interval_nt operator/(const CGAL_Interval_nt& d) const
  {
    CGAL_Interval_nt tmp;
    CGAL_FPU_set_rounding_to_infinity();
    tmp = ((CGAL_Interval_nt_advanced) *this) / (CGAL_Interval_nt_advanced) d;
    CGAL_FPU_set_rounding_to_nearest();
    return tmp;
  }
};


inline CGAL_Interval_nt_advanced sqrt(const CGAL_Interval_nt_advanced& d)
{
  CGAL_Interval_nt_advanced tmp;
  CGAL_FPU_set_rounding_to_minus_infinity();
  tmp.inf = -sqrt(-d.inf);
  CGAL_FPU_set_rounding_to_infinity();
  tmp.sup = sqrt(d.sup);
  return tmp;
}

inline CGAL_Interval_nt sqrt(const CGAL_Interval_nt& d)
{
  CGAL_Interval_nt tmp = sqrt( (CGAL_Interval_nt_advanced) d);
  CGAL_FPU_set_rounding_to_nearest();
  return tmp;
}

inline double CGAL_to_double(const CGAL_Interval_nt_advanced& d)
{
    return (d.sup-d.inf)*.5;
}

ostream& operator<<(ostream& os, CGAL_Interval_nt_advanced& d)
{
  os << "[" << d.lower_bound() << ";" << d.upper_bound() << "]";
  return os;
}

#endif /* CGAL_INTERVAL_ARITHMETIC_H */
