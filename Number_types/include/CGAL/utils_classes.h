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

#ifdef CGAL_USE_SSE2_MAX
#include <CGAL/sse2.h>
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

#ifdef CGAL_USE_SSE2_MAX

  inline double sse2max(double a, double b, double c, double d)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);
  __m128d C =_mm_load_sd(&c);
  __m128d D =_mm_load_sd(&d);

  __m128d AB = _mm_max_sd(A,B);
  __m128d CD = _mm_max_sd(C,D);
  A = _mm_max_sd(AB,CD);
  _mm_store_sd(&a, A);
  return a;
}

inline double sse2max(double a, double b, double c)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);
  __m128d C =_mm_load_sd(&c);

  __m128d AB = _mm_max_sd(A,B);
  A = _mm_max_sd(AB,C);
  _mm_store_sd(&a, A);
  return a;
}

inline double sse2max(double a, double b)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);

  __m128d C = _mm_max_sd(A,B);
  _mm_store_sd(&a, C);
  return a;
}


#if 0
// Doing things in parallel seems the way to go
// but copying to/from arrays has too much overhead
//  a = max(a,a2) b = max(b,b2)
inline void sse2mmax2(double& a, double a2, double& b, double b2)
{
   CGAL_ALIGN_16 double res[2];
  res[0] = a;
  res[1] = b;
  __m128d F =_mm_load_pd(res);
  res[0] = a2;
  res[1] = b2;
  __m128d S =_mm_load_pd(res);

  __m128d C = _mm_max_pd(F,S);

  _mm_store_pd(res, C);
  a = res[0];
  b = res[1];
}
#endif


  inline double sse2min(double a, double b, double c, double d)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);
  __m128d C =_mm_load_sd(&c);
  __m128d D =_mm_load_sd(&d);

  __m128d AB = _mm_min_sd(A,B);
  __m128d CD = _mm_min_sd(C,D);
  A = _mm_min_sd(AB,CD);
  _mm_store_sd(&a, A);
  return a;
}

inline double sse2min(double a, double b, double c)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);
  __m128d C =_mm_load_sd(&c);

  __m128d AB = _mm_min_sd(A,B);
  A = _mm_min_sd(AB,C);
  _mm_store_sd(&a, A);
  return a;
}

inline double sse2min(double a, double b)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);

  __m128d C = _mm_min_sd(A,B);
  _mm_store_sd(&a, C);
  return a;
}

inline void sse2minmax(double& a, double b, double& c)
{
  __m128d A =_mm_load_sd(&a);
  __m128d B =_mm_load_sd(&b);
  __m128d C =_mm_load_sd(&c);

  __m128d AB = _mm_min_sd(A,B);
  A = _mm_min_sd(AB,C);
  _mm_store_sd(&a, A);

  AB = _mm_max_pd(A,B);
  C = _mm_max_sd(AB,C);
  _mm_store_sd(&c, C);
}

#endif // CGAL_USE_SSE2_MAX

template <>
struct Max<double> :public CGAL::cpp98::binary_function< double, double, double > {
 Max() {}

 double operator()( const double& x, const double& y) const
    {
#ifdef CGAL_USE_SSE2_MAX
      return sse2max(x,y);
#else
      return (std::max)( x, y);
#endif
 }

  double operator()( double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return sse2max(x,y,z);
#else
    return (std::max)((std::max)( x, y), z);
#endif
  }

  double operator()( double w,double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return sse2max(w,x,y,z);
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
      return sse2min(x,y);
#else
      return (std::min)( x, y);
#endif
 }

  double operator()( double x, double y, double z) const
  {
#ifdef CGAL_USE_SSE2_MAX
    return sse2min(x,y,z);
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

} //namespace CGAL

#endif // CGAL_UTILS_CLASSES_H
