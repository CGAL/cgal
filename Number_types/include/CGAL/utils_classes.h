// Copyright (c) 2006-2008  Max-Planck-Institute Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-sb.mpg.de>

#ifndef CGAL_UTILS_CLASSES_H
#define CGAL_UTILS_CLASSES_H

#include <CGAL/config.h>
#include <functional> // for std::less
#include <algorithm>  // for std::min and max

#if defined(_MSC_VER) && (defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#include <emmintrin.h>
#define CGAL_USE_SSE2_MINMAX
#endif

namespace CGAL {

template < class A, class B = A >
struct Equal_to : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x == y; }
};

template < class A, class B = A >
struct Not_equal_to : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x != y; }
};

template < class A, class B = A >
struct Greater : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x > y; }
};

template < class A, class B = A >
struct Less : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x < y; }
};

template < class A, class B = A >
struct Greater_equal : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x >= y; }
};

template < class A, class B = A >
struct Less_equal : public CGAL::cpp98::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x <= y; }
};

template < class NT, class Less = std::less< NT > >
struct Min :public CGAL::cpp98::binary_function< NT, NT, NT > {
 Min() {}
 Min(const Less& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::min)( x, y, c); }
protected:
 Less c;
};

template < class NT, class Less = std::less< NT > >
struct Max :public CGAL::cpp98::binary_function< NT, NT, NT > {
 Max() {}
 Max(const Less& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::max)( x, y, c); }
protected:
 Less c;
};

// Both return exactly std::max(a, b) / std::min(a, b), NaN behaviour included.
inline double dmax(double a, double b)
{
#ifdef CGAL_USE_SSE2_MINMAX
  return _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd(b), _mm_set_sd(a)));   // (a < b) ? b : a
#else
  return (std::max)(a, b);
#endif
}

inline double dmin(double a, double b)
{
#ifdef CGAL_USE_SSE2_MINMAX
  return _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd(b), _mm_set_sd(a)));   // (b < a) ? b : a
#else
  return (std::min)(a, b);
#endif
}

inline double dmax(double a, double b, double c)
{
  return dmax(dmax(a, b), c);
}

inline double dmax(double a, double b, double c, double d)
{
  return dmax(dmax(a, b), dmax(c, d));   // tree: depth 2 instead of 3
}

inline double dmax(double a, double b, double c, double d, double e)
{
  return dmax(dmax(a, b, c, d), e);   // depth 3, the minimum for 5 values
}

inline double dmin(double a, double b, double c)
{
  return dmin(dmin(a, b), c);
}

inline double dmin(double a, double b, double c, double d)
{
  return dmin(dmin(a, b), dmin(c, d));
}

inline double dmin(double a, double b, double c, double d, double e)
{
  return dmin(dmin(a, b, c, d), e);
}

template <>
struct Max<double> :public CGAL::cpp98::binary_function< double, double, double > {
 Max() {}

 double operator()( const double& x, const double& y) const
    {
#ifdef CGAL_USE_SSE2_MAX
      return dmax(x,y);
#else
      return (std::max)( x, y);
#endif
 }

  double operator()( double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return dmax(x,y,z);
#else
    return (std::max)((std::max)( x, y), z);
#endif
  }

  double operator()( double w,double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return dmax(w,x,y,z);
#else
    return (std::max)((std::max)( x, y), (std::max)(w,z));
#endif
  }
};

template <>
struct Min<double> :public CGAL::cpp98::binary_function< double, double, double > {
 Min() {}

 double operator()( const double& x, const double& y) const
    {
#ifdef CGAL_USE_SSE2_MAX
      return dmin(x,y);
#else
      return (std::min)( x, y);
#endif
 }

  double operator()( double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return dmin(x,y,z);
#else
    return (std::min)((std::min)( x, y), z);
#endif
  }

  double operator()( double w,double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return sse2min(w,x,y,z);
#else
    return (std::min)((std::min)( x, y), (std::min)(w,z));
#endif
  }
};
template< class T >
class Is_valid
  : public CGAL::cpp98::unary_function< T, bool > {
  public:
    bool operator()( const T& ) const {
      return true;
    };
};

namespace internal
{
// utility class to be used for calling exact(Lazy) when doing accumulation with EPECK
template <class NT>
struct Evaluate
{
  template <class T>
  void operator()(const T&)
  {}
};
} // internal namespace

} //namespace CGAL

#endif // CGAL_UTILS_CLASSES_H
