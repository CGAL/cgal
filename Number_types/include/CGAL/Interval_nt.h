// Copyright (c) 1998-2019
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion, Michael Hemmer, Marc Glisse

#ifndef CGAL_INTERVAL_NT_H
#define CGAL_INTERVAL_NT_H

// This file contains the description of the following classes:
// - Interval_nt<false>  It's a number type that needs the FPU rounding mode
//                       to be set to +inf.  It is also typedef'd to
//                       Interval_nt_advanced for backward compatibility.
// - Interval_nt<true>   Same but it does the rounding mode itself so you
//                       don't have to worry about it.  But it's slower.
//
// Note: When rounding is towards +infinity, to make an operation rounded
// towards -infinity, it's enough to take the opposite of some of the operand,
// and the opposite of the result (see operator+, operator*,...).

// TODO :
// - test whether stopping constant propagation only in functions taking
//   double as arguments, improves performance.

#include <utility> // for std::pair
#include <CGAL/number_type_config.h>
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/number_utils.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/double.h>
#include <CGAL/FPU.h>
#include <CGAL/IO/io.h>
#include <iostream>
#include <boost/operators.hpp>

#ifdef __GNUC__
// gcc's __builtin_constant_p does not like arguments with side effects. Be
// careful not to use this macro for something that the compiler will have
// trouble eliminating as dead code.
# define CGAL_CST_TRUE(X) ({ bool _ugly_ = (X); __builtin_constant_p(_ugly_) && _ugly_; })
#else
# define CGAL_CST_TRUE(X) false
#endif

namespace CGAL {

template <bool Protected = true>
class Interval_nt
{
  typedef Interval_nt<Protected>     IA;
  typedef std::pair<double, double>  Pair;

public:

  typedef double      value_type;

  typedef Uncertain_conversion_exception            unsafe_comparison;
  typedef Checked_protect_FPU_rounding<Protected>   Internal_protector;
  typedef Protect_FPU_rounding<!Protected>          Protector;

  Interval_nt()
#ifndef CGAL_NO_ASSERTIONS
# ifdef CGAL_USE_SSE2
      : val(_mm_setr_pd(-1, 0))
# else
      : _inf(-1), _sup(0)
# endif
             // to early and deterministically detect use of uninitialized
#endif
    {}

  Interval_nt(int i)
  { *this = static_cast<double>(i); }

  Interval_nt(unsigned i)
  { *this = static_cast<double>(i); }

  Interval_nt(long long i)
  {
    // gcc ignores -frounding-math when converting integers to floats.
    // Is this safe against excess precision? -- Marc Glisse, Dec 2012
    double d = static_cast<double>(i);
    *this = d;
#ifdef __GNUC__
    long long safe = 1LL << 52; // Use numeric_limits?
    bool exact = ((long long)d == i) || (i <= safe && i >= -safe);
    if (!CGAL_CST_TRUE(exact))
#endif
      *this += smallest();
  }

  Interval_nt(unsigned long long i)
  {
    double d = static_cast<double>(i);
    *this = d;
#ifdef __GNUC__
    unsigned long long safe = 1ULL << 52; // Use numeric_limits?
    bool exact = ((unsigned long long)d == i) || (i <= safe);
    if (!CGAL_CST_TRUE(exact))
#endif
      *this += smallest();
  }

  Interval_nt(long i)
  {
    *this = (sizeof(int)==sizeof(long)) ?
      Interval_nt((int)i) :
      Interval_nt((long long)i);
  }

  Interval_nt(unsigned long i)
  {
    *this = (sizeof(int)==sizeof(long)) ?
      Interval_nt((unsigned)i) :
      Interval_nt((unsigned long long)i);
  }

  Interval_nt(double d)
  {
    CGAL_assertion(is_finite(d));
    *this = Interval_nt(d, d);
  }

// The Intel compiler on Linux is aggressive with constant propagation and
// it seems there is no flag to stop it, so disable this check for it.
#if !defined(CGAL_DISABLE_ROUNDING_MATH_CHECK) && \
    defined(__INTEL_COMPILER) && defined(__linux)
#  define CGAL_DISABLE_ROUNDING_MATH_CHECK
#endif

#ifdef CGAL_USE_SSE2
  // This constructor should really be private, like the simd() function, but
  // that would mean a lot of new friends, so they are only undocumented.
  explicit Interval_nt(__m128d v) : val(v) {}
#endif

  Interval_nt(double i, double s)
#ifdef CGAL_USE_SSE2
    : val(_mm_setr_pd(-i, s))
#else
    : _inf(-i), _sup(s)
#endif
  {
    // Previously it was:
    //    CGAL_assertion_msg(!(i>s);
    // But MSVC++ 2012 optimizes the test "!(i>s)" to "i<=s", even with
    // /fp:strict. If 'i' or 's' is a NaN, that makes a difference.
    CGAL_assertion_msg( (!is_valid(i)) || (!is_valid(s)) || (!(i>s)),
              " Variable used before being initialized (or CGAL bug)");
#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
    CGAL_assertion_code((void) tester;) // Necessary to trigger a runtime test of rounding modes.
#endif
  }

  Interval_nt(const Pair & p)
  { *this = Interval_nt(p.first, p.second); }

  IA operator-() const
  {
#ifdef CGAL_USE_SSE2
    return IA (swap_m128d(val));
#else
    return IA (-sup(), -inf());
#endif
  }

  IA & operator+= (const IA &d) { return *this = *this + d; }
  IA & operator-= (const IA &d) { return *this = *this - d; }
  IA & operator*= (const IA &d) { return *this = *this * d; }
  IA & operator/= (const IA &d) { return *this = *this / d; }

  bool is_point() const
  {
    return sup() == inf();
  }

  bool is_same (const IA & d) const
  {
#ifdef CGAL_USE_SSE2
    // Faster to answer yes, but slower to answer no.
    return _mm_movemask_pd (_mm_cmpneq_pd (val, d.val)) == 0;
#else
    return inf() == d.inf() && sup() == d.sup();
#endif
  }

  bool do_overlap (const IA & d) const
  {
#ifdef CGAL_USE_SSE2
    __m128d m = _mm_set1_pd (-0.);
    __m128d y = _mm_xor_pd ((-d).val, m); // {-ds,di}
    __m128d c = _mm_cmplt_pd (val, y); // {i>ds,s<di}
    return _mm_movemask_pd (c) == 0;
#else
    return !(d.inf() > sup() || d.sup() < inf());
#endif
  }

  double inf() const
  {
#ifdef CGAL_USE_SSE2
    return -_mm_cvtsd_f64(val);
#else
    return -_inf;
#endif
  }
  double sup() const
  {
#ifdef CGAL_USE_SSE2
    return _mm_cvtsd_f64(swap_m128d(val));
    // The following is a bit more natural, but
    // - it is too opaque
    // - it is a less likely CSE candidate
    // return _mm_cvtsd_f64(_mm_unpackhi_pd(val, val));
#else
    return _sup;
#endif
  }
#ifdef CGAL_USE_SSE2
  __m128d simd() const { return val; }
#endif

  std::pair<double, double> pair() const
  {
    return std::pair<double, double>(inf(), sup());
  }

  static IA largest()
  {
    return IA(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  }

  static IA smallest()
  {
    return IA(-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE);
  }

#if 0 // def CGAL_HISTOGRAM_PROFILER  // not yet ready
  ~Interval_nt()
  {
    CGAL_HISTOGRAM_PROFILER("[Interval_nt relative precision in log2 scale]",
                             (unsigned) ( ::log(relative_precision(*this))) / ::log(2.0) )  );
  }
#endif

private:
  // Pair inf_sup;
  // The value stored in _inf is the negated lower bound.
  // TODO: experiment with different orders of the values in the SSE2 register,
  // for instance {sup, -inf}, or {inf, -sup}, and adapt users to query the low
  // value in priority. {-inf, sup} has the drawback that neither inf nor sup
  // is free to access.
#ifdef CGAL_USE_SSE2
  __m128d val;
#else
  double _inf, _sup;
#endif

  struct Test_runtime_rounding_modes {
    Test_runtime_rounding_modes()
    {
      // We test whether GCC's -frounding-math option has been forgotten.
      // The macros CGAL_IA_MUL and CGAL_IA_DIV stop constant propagation only
      // on the second argument, so if -fno-rounding-math, the compiler optimizes
      // the 2 negations and we get wrong rounding.
      typename Interval_nt<>::Internal_protector P;
      CGAL_assertion_msg(-CGAL_IA_MUL(-1.1, 10.1) != CGAL_IA_MUL(1.1, 10.1),
                         "Wrong rounding: did you forget the  -frounding-math  option if you use GCC (or  -fp-model strict  for Intel)?");
      CGAL_assertion_msg(-CGAL_IA_DIV(-1., 10) != CGAL_IA_DIV(1., 10),
                         "Wrong rounding: did you forget the  -frounding-math  option if you use GCC (or  -fp-model strict  for Intel)?");
    }
  };

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  static const Test_runtime_rounding_modes tester;
#endif

  friend
    Uncertain<bool>
    operator<(const Interval_nt &a, const Interval_nt &b)
    {
      if (a.sup()  < b.inf()) return true;
      if (a.inf() >= b.sup()) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>(const Interval_nt &a, const Interval_nt &b)
    { return b < a; }

  friend
    Uncertain<bool>
    operator<=(const Interval_nt &a, const Interval_nt &b)
    {
      if (a.sup() <= b.inf()) return true;
      if (a.inf() >  b.sup()) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>=(const Interval_nt &a, const Interval_nt &b)
    { return b <= a; }

  friend
    Uncertain<bool>
    operator==(const Interval_nt &a, const Interval_nt &b)
    {
      if (b.inf() >  a.sup() || b.sup() <  a.inf()) return false;
      if (b.inf() == a.sup() && b.sup() == a.inf()) return true;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator!=(const Interval_nt &a, const Interval_nt &b)
    { return ! (a == b); }


  // Mixed operators with double.

  friend
    Uncertain<bool>
    operator<(double a, const Interval_nt &b)
    {
      if (a < b.inf()) return true;
      if (a >= b.sup()) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>(double a, const Interval_nt &b)
    { return b < a; }

  friend
    Uncertain<bool>
    operator<=(double a, const Interval_nt &b)
    {
      if (a <= b.inf()) return true;
      if (a >  b.sup()) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>=(double a, const Interval_nt &b)
    { return b <= a; }

  friend
    Uncertain<bool>
    operator==(double a, const Interval_nt &b)
    {
      if (b.inf() >  a || b.sup() <  a) return false;
      if (b.is_point()) return true;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator!=(double a, const Interval_nt &b)
    { return ! (a == b); }

  friend
    Uncertain<bool>
    operator<(const Interval_nt &a, double b)
    {
      if (a.sup() <  b) return true;
      if (a.inf() >= b) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>(const Interval_nt &a, double b)
    { return b < a; }

  friend
    Uncertain<bool>
    operator<=(const Interval_nt &a, double b)
    {
      if (a.sup() <= b) return true;
      if (a.inf() >  b) return false;
      return Uncertain<bool>::indeterminate();
    }

  friend
    Uncertain<bool>
    operator>=(const Interval_nt &a, double b)
    { return b <= a; }

  friend
    Uncertain<bool>
    operator==(const Interval_nt &a, double b)
    {
      return b == a;
    }

  friend
    Uncertain<bool>
    operator!=(const Interval_nt &a, double b)
    { return b != a; }

  friend
    std::ostream & operator<< (std::ostream &os, const Interval_nt & I )
    {
      return os << "[" << I.inf() << ";" << I.sup() << "]";
    }

#define CGAL_SWALLOW(IS,CHAR)                    \
  {                                              \
    char c;                                      \
    do is.get(c); while (isspace(c));            \
    if (c != CHAR) {                             \
      is.setstate(std::ios_base::failbit);       \
    }                                            \
  }                                              \

  friend
    std::istream & operator>> (std::istream &is, Interval_nt & I)
    {
      char c;
      do is.get(c); while (isspace(c));
      is.putback(c);
      if(c == '['){ // read original output from operator <<
        double inf,sup;
        CGAL_SWALLOW(is, '[');// read the "["
        is >> iformat(inf);
        CGAL_SWALLOW(is, ';');// read the ";"
        is >> iformat(sup);
        CGAL_SWALLOW(is, ']');// read the "]"
        I = Interval_nt(inf,sup);
      }else{ //read double (backward compatibility)
        double d;
        is >> d;
        I = d;
      }
      return is;
    }
#undef CGAL_SWALLOW

  friend
    Interval_nt
    operator+ (const Interval_nt &a, const Interval_nt & b)
    {
      Internal_protector P;
#ifdef CGAL_USE_SSE2
      __m128d aa = IA_opacify128(a.simd());
      __m128d bb = IA_opacify128_weak(b.simd());
      __m128d r = _mm_add_pd(aa, bb);
      return Interval_nt(IA_opacify128(r));
#else
      return Interval_nt (-CGAL_IA_ADD(-a.inf(), -b.inf()),
          CGAL_IA_ADD(a.sup(), b.sup()));
#endif
    }

  // MSVC does not define __SSE3__
#if defined CGAL_USE_SSE2 && (defined __SSE3__ || defined __AVX__)
  friend
    Interval_nt
    operator+ (double a, const Interval_nt & b)
    {
      Internal_protector P;
      __m128d aa = _mm_set1_pd(IA_opacify(a));
      __m128d bb = IA_opacify128_weak(b.simd());
      __m128d r = _mm_addsub_pd(bb, aa);
      return Interval_nt(IA_opacify128(r));
    }

  friend
    Interval_nt
    operator+ (const Interval_nt & a, double b)
    {
      return b + a;
    }
#endif

  friend
    Interval_nt
    operator+( const Interval_nt& a ) {
      return a;
    }

  friend
    Interval_nt
    operator- (const Interval_nt &a, const Interval_nt & b)
    {
#ifdef CGAL_USE_SSE2
      return a+-b;
#else
      Internal_protector P;
      return Interval_nt(-CGAL_IA_ADD(b.sup(), -a.inf()),
          CGAL_IA_ADD(a.sup(), -b.inf()));
#endif
    }

#ifdef CGAL_USE_SSE2
  friend
    Interval_nt
    operator- (double a, const Interval_nt & b)
    {
      return a+-b;
    }
#endif

  friend
    Interval_nt
    operator* (const Interval_nt &a, const Interval_nt & b)
    {
#if 0
      // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=88626
      if(CGAL_CST_TRUE(a.is_point()))
        return a.inf() * b;
      else if(CGAL_CST_TRUE(b.is_point()))
        return a * b.inf();
#endif
      Internal_protector P;
#ifdef CGAL_USE_SSE2
# if !defined __SSE4_1__ && !defined __AVX__
      // Brutal, compute all products in all directions.
      // The actual winner (by a hair) on recent hardware before removing NaNs.
      __m128d aa = IA_opacify128_weak(a.simd());            // {-ai,as}
    __m128d bb = b.simd();                                // {-bi,bs}
    __m128d m = _mm_set_sd(-0.);                          // {-0,+0}
    __m128d m1 = _mm_set1_pd(-0.);                        // {-0,-0}
    __m128d ax = swap_m128d (aa);                         // {as,-ai}
    __m128d ap = _mm_xor_pd (ax, m1);                     // {-as,ai}
    __m128d bz = _mm_xor_pd(bb, m);                       // {bi,bs}
    bz = IA_opacify128(bz);
    __m128d c = swap_m128d (bz);                          // {bs,bi}

    // The multiplications could produce some NaN, with 0 * inf. Replacing it with inf is safe.
    // min(x,y) (the order is essential) returns its second argument when the first is NaN.
    // An IEEE 754-2019 maximum could help.
    __m128d big = IA::largest().simd();
    __m128d x1 = _mm_mul_pd(aa,bz);                       // {-ai*bi,as*bs}
    //x1 = _mm_min_pd(x1,big); // no NaN
    __m128d x2 = _mm_mul_pd(aa,c);                        // {-ai*bs,as*bi}
    x2 = _mm_min_pd(x2,big); // no NaN
    __m128d x3 = _mm_mul_pd(ap,bz);                       // {-as*bi,ai*bs}
    //x3 = _mm_min_pd(x3,big); // no NaN
    __m128d x4 = _mm_mul_pd(ap,c);                        // {-as*bs,ai*bi}
    x4 = _mm_min_pd(x4,big); // no NaN

    __m128d y1 = _mm_max_pd(x1,x2);
    __m128d y2 = _mm_max_pd(x3,x4);
    __m128d r = _mm_max_pd (y1, y2);
    // Alternative with fewer instructions but more dependency
    // __m128d r = _mm_max_pd(x1,_mm_max_pd(x2,_mm_max_pd(x3,_mm_min_pd(x4,big))));
    return IA (IA_opacify128(r));
# elif 1
    // we want to multiply ai,as with {ai<0?-bs:-bi,as<0?bi:bs}
    // we want to multiply as,ai with {as<0?-bs:-bi,ai<0?bi:bs}
    // requires SSE4 (otherwise use _mm_cmplt_pd, _mm_and_pd, _mm_andnot_pd and _mm_or_pd to avoid blendv)
    // probably faster on older hardware
    __m128d m = _mm_set_sd(-0.);                          // {-0,+0}
    __m128d m1 = _mm_set1_pd(-0.);                        // {-0,-0}
    __m128d big = IA::largest().simd();
    __m128d aa = a.simd();                                // {-ai,as}
    __m128d az = _mm_xor_pd(aa, m);                       // {ai,as}
    az = IA_opacify128_weak(az);
    __m128d azp = swap_m128d (az);                        // {as,ai}
    __m128d bb = IA_opacify128(b.simd());                 // {-bi,bs}
    __m128d bx = swap_m128d (bb);                         // {bs,-bi}
    __m128d bp = _mm_xor_pd(bx, m1);                      // {-bs,bi}
    __m128d x = _mm_blendv_pd (bb, bp, az);               // {ai<0?-bs:-bi,as<0?bi:bs}
    __m128d y = _mm_blendv_pd (bb, bp, azp);              // {as<0?-bs:-bi,ai<0?bi:bs}
    __m128d p1 = _mm_mul_pd (az, x);
    //p1 = _mm_min_pd(p1,big); // no NaN
    __m128d p2 = _mm_mul_pd (azp, y);
    p2 = _mm_min_pd(p2,big); // no NaN
    __m128d r = _mm_max_pd (p1, p2);
    return IA (IA_opacify128(r));
# elif 0
    // we want to multiply -ai,as with {ai>0?bi:bs,as<0?bi:bs}
    // we want to multiply -as,ai with {as<0?bs:bi,ai>0?bs:bi}
    // slightly worse than the previous one
    __m128d m1 = _mm_set1_pd(-0.);                        // {-0,-0}
    __m128d big = IA::largest().simd();
    __m128d aa = IA_opacify128_weak(a.simd());            // {-ai,as}
    __m128d ax = swap_m128d (aa);                         // {as,-ai}
    __m128d ap = _mm_xor_pd (ax, m1);                     // {-as,ai}
    __m128d bb = IA_opacify128(b.simd());                 // {-bi,bs}
    double bi = -_mm_cvtsd_f64(bb);
    double bs = _mm_cvtsd_f64(_mm_unpackhi_pd(bb,bb));
    __m128d bbi = _mm_set1_pd(bi);                        // {bi,bi}
    __m128d bbs = _mm_set1_pd(bs);                        // {bs,bs}
    __m128d x = _mm_blendv_pd (bbs, bbi, aa);             // {ai>0?bi:bs,as<0?bi:bs}
    __m128d y = _mm_blendv_pd (bbi, bbs, ax);             // {as<0?bs:bi,ai>0?bs:bi}
    __m128d p1 = _mm_mul_pd (aa, x);
    //p1 = _mm_min_pd(p1,big); // no NaN
    __m128d p2 = _mm_mul_pd (ap, y);
    p2 = _mm_min_pd(p2,big); // no NaN
    __m128d r = _mm_max_pd (p1, p2);
    return IA (IA_opacify128(r));
# else
    // AVX version of the brutal method, same running time or slower
    __m128d aa = IA_opacify128_weak(a.simd());            // {-ai,as}
    __m128d bb = b.simd();                                // {-bi,bs}
    __m256d big = _mm256_set1_pd(std::numeric_limits<double>::infinity());
    __m128d m = _mm_set_sd(-0.);                          // {-0,+0}
    __m128d m1 = _mm_set1_pd(-0.);                        // {-0,-0}
    __m128d ax = swap_m128d (aa);                         // {as,-ai}
    __m128d ap = _mm_xor_pd (ax, m1);                     // {-as,ai}
    __m128d bz = _mm_xor_pd(bb, m);                       // {bi,bs}
    bz = IA_opacify128(bz);
    __m256d X = _mm256_set_m128d(ap,aa);                  // {-ai,as,-as,ai}
    __m256d Y1 = _mm256_set_m128d(bz,bz);                 // {bi,bs,bi,bs}
    __m256d Y2 = _mm256_permute_pd(Y1,5);                 // {bs,bi,bs,bi}
    __m256d Z1 = _mm256_mul_pd(X,Y1);
    //Z1 = _mm256_min_pd(Z1,big); // no NaN
    __m256d Z2 = _mm256_mul_pd(X,Y2);
    Z2 = _mm256_min_pd(Z2,big); // no NaN
    __m256d Z = _mm256_max_pd(Z1,Z2);
    __m128d z1 = _mm256_castpd256_pd128(Z);
    __m128d z2 = _mm256_extractf128_pd(Z,1);
    __m128d r = _mm_max_pd (z1, z2);
    return IA (IA_opacify128(r));
# endif
#else
    // TODO: try to move some NaN tests out of the hot path (test a.inf()>0 instead of >=0?).
    if (a.inf() >= 0.0)                                        // a>=0
    {
      // b>=0     [a.inf()*b.inf(); a.sup()*b.sup()]
      // b<=0     [a.sup()*b.inf(); a.inf()*b.sup()]
      // b~=0     [a.sup()*b.inf(); a.sup()*b.sup()]
      double aa = a.inf(), bb = a.sup();
      if (bb <= 0.) return 0.; // In case b has an infinite bound, avoid NaN.
      if (b.inf() < 0.0)
      {
        aa = bb;
        if (b.sup() < 0.0)
          bb = a.inf();
      }
      double r = (b.sup() == 0) ? 0. : CGAL_IA_MUL(bb, b.sup()); // In case bb is infinite, avoid NaN.
      return IA(-CGAL_IA_MUL(aa, -b.inf()), r);
    }
    else if (a.sup()<=0.0)                                // a<=0
    {
      // b>=0     [a.inf()*b.sup(); a.sup()*b.inf()]
      // b<=0     [a.sup()*b.sup(); a.inf()*b.inf()]
      // b~=0     [a.inf()*b.sup(); a.inf()*b.inf()]
      double aa = a.sup(), bb = a.inf();
      if (b.inf() < 0.0)
      {
        aa=bb;
        if (b.sup() <= 0.0)
          bb=a.sup();
      }
      else if (b.sup() <= 0) return 0.; // In case a has an infinite bound, avoid NaN.
      return IA(-CGAL_IA_MUL(-bb, b.sup()), CGAL_IA_MUL(-aa, -b.inf()));
    }
    else                                                // 0 \in a
    {
      if (b.inf()>=0.0) {                                // b>=0
        if (b.sup()<=0.0)
          return 0.; // In case a has an infinite bound, avoid NaN.
        else
          return IA(-CGAL_IA_MUL(-a.inf(), b.sup()),
              CGAL_IA_MUL( a.sup(), b.sup()));
      }
      if (b.sup()<=0.0) {                                // b<=0
        return IA(-CGAL_IA_MUL( a.sup(), -b.inf()),
            CGAL_IA_MUL(-a.inf(), -b.inf()));
      }
      // 0 \in b
      double tmp1 = CGAL_IA_MUL(-a.inf(),  b.sup());
      double tmp2 = CGAL_IA_MUL( a.sup(), -b.inf());
      double tmp3 = CGAL_IA_MUL(-a.inf(), -b.inf());
      double tmp4 = CGAL_IA_MUL( a.sup(),  b.sup());
      return IA(-(std::max)(tmp1,tmp2), (std::max)(tmp3,tmp4));
    }
#endif
    }

  friend
    Interval_nt
    operator* (double a, Interval_nt b)
    {
      CGAL_assertion(is_finite(a));
      // return Interval_nt(a)*b;
      Internal_protector P;
      if (a < 0) { a = -a; b = -b; }
      // Now a >= 0
#ifdef CGAL_USE_SSE2
      // TODO: try/benchmark a branchless version
      __m128d bb = IA_opacify128_weak(b.simd());
      __m128d aa = _mm_set1_pd(IA_opacify(a));
      __m128d r = _mm_mul_pd(aa, bb);
      // In case a is 0 and b has an infinite bound. This returns an interval
      // larger than necessary, but is likely faster to produce.
      r = _mm_min_pd(r,largest().simd());
      return IA(IA_opacify128(r));
#else
      else if (!(a > 0)) return 0.; // We could test this before the SSE block and remove the minpd line.
      return IA(-CGAL_IA_MUL(a, -b.inf()), CGAL_IA_MUL(a, b.sup()));
#endif
    }

  friend
    Interval_nt
    operator* (const Interval_nt & a, double b)
    {
      return b * a;
    }

  friend
    Interval_nt
    operator/ (const Interval_nt &a, const Interval_nt & b)
    {
#if 0
      // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=88626
      if(CGAL_CST_TRUE(a.is_point()))
        return a.inf() / b;
      else if(CGAL_CST_TRUE(b.is_point()))
        return a / b.inf();
#endif
      Internal_protector P;
#if defined CGAL_USE_SSE2 && (defined __SSE4_1__ || defined __AVX__)
      //// not a tight bound, but easy:
      // return CGAL::inverse(b)*a;
# if 1
      // Current fastest
      // if b>0 we want [ai/(ai>0?bs:bi),as/(as>0?bi:bs)]
      // if b<0 we want [as/(as>0?bs:bi),ai/(ai>0?bi:bs)]
      __m128d m = _mm_set_sd(-0.);
      __m128d aa = a.simd();
      __m128d bb = b.simd();
      int i = _mm_movemask_pd(_mm_cmpge_pd(bb, _mm_set1_pd(0.)));
      if(i==3) return largest(); // bi<=0 && bs>=0
      __m128d ap = _mm_xor_pd(aa, m); // {ai, as}
    __m128d ax = swap_m128d(ap); // {as, ai}
    __m128d bp = _mm_xor_pd(bb, m); // {bi, bs}
    __m128d bx = swap_m128d(bp); // {bs, bi}
    __m128d num = _mm_blendv_pd(ap, ax, bp); // {(b>0)?ai:as, (b>0)?as:ai}
    __m128d d = _mm_blendv_pd(bx, bp, num);
    // Can we rearrange things so we need fewer xor?
    __m128d den = _mm_xor_pd(d, m);
    num = IA_opacify128_weak(num);
    den = IA_opacify128(den);
    __m128d r = _mm_div_pd(num, den);
    return IA (IA_opacify128(r));
# else
    // Similar to the multiplication, but slow, because divisions are slow
    // if b>0 we want [-max(-ai/bi,-ai/bs),max(as/bi,as/bs)] {-ai,as}/{bi,bs} {-ai,as}/{bs,bi}
    // if b<0 we want [-max(-as/bi,-as/bs),max(ai/bi,ai/bs)] {-as,ai}/{bi,bs} {-as,ai}/{bs,bi}
    __m128d m = _mm_set_sd(-0.);
    __m128d m1 = _mm_set1_pd(-0.);
    __m128d aa = a.simd(); // {-ai, as}
    __m128d bb = b.simd(); // {-bi, bs}
    int i = _mm_movemask_pd(_mm_cmpge_pd(bb, _mm_set1_pd(0.)));
    if(i==3) return largest(); // bi<=0 && bs>=0
    __m128d ap = _mm_xor_pd(aa, m1); // {ai, -as}
    __m128d ax = swap_m128d(ap); // {-as, ai}
    __m128d bp = _mm_xor_pd(bb, m); // {bi, bs}
    __m128d num = _mm_blendv_pd(aa, ax, bp);
    num = IA_opacify128_weak(num);
    bp = IA_opacify128(bp);
    __m128d bx = swap_m128d(bp); // {bs, bi}
    __m128d d1 = _mm_div_pd(num, bp);
    __m128d d2 = _mm_div_pd(num, bx);
    __m128d r = _mm_max_pd(d1, d2);
    return IA (IA_opacify128(r));
# endif
#else
    if (b.inf() > 0.0)                                // b>0
    {
      // e>=0        [a.inf()/b.sup(); a.sup()/b.inf()]
      // e<=0        [a.inf()/b.inf(); a.sup()/b.sup()]
      // e~=0        [a.inf()/b.inf(); a.sup()/b.inf()]
      double aa = b.sup(), bb = b.inf();
      if (a.inf() < 0.0)
      {
        aa = bb;
        if (a.sup() < 0.0)
          bb = b.sup();
      }
      return IA(-CGAL_IA_DIV(-a.inf(), aa), CGAL_IA_DIV(a.sup(), bb));
    }
    else if (b.sup()<0.0)                        // b<0
    {
      // e>=0        [a.sup()/b.sup(); a.inf()/b.inf()]
      // e<=0        [a.sup()/b.inf(); a.inf()/b.sup()]
      // e~=0        [a.sup()/b.sup(); a.inf()/b.sup()]
      double aa = b.sup(), bb = b.inf();
      if (a.inf() < 0.0)
      {
        bb = aa;
        if (a.sup() < 0.0)
          aa = b.inf();
      }
      return IA(-CGAL_IA_DIV(a.sup(), -aa), CGAL_IA_DIV(a.inf(), bb));
    }
    else                                        // b~0
      return largest();
    // We could do slightly better -> [0;infinity] when b.sup()==0,
    // but is this worth ?
#endif
    }

  // Without SSE2, let it use the function above.
#ifdef CGAL_USE_SSE2
  friend
    Interval_nt
    operator/ (double a, const Interval_nt & b)
    {
      int i = _mm_movemask_pd(_mm_cmpge_pd(b.simd(), _mm_set1_pd(0.)));
      if(i==3) return largest(); // bi<=0 && bs>=0
      __m128d aa, xx;
      if(a>0){
        aa = _mm_set1_pd(-a);
        xx = (-b).simd();
      } else if(a<0){
        aa = _mm_set1_pd(a);
        xx = b.simd();
      } else return 0.;
      Internal_protector P;
      __m128d r = _mm_div_pd(IA_opacify128_weak(aa), IA_opacify128(xx));
      return Interval_nt(IA_opacify128(r));
    }

  friend
    Interval_nt
    operator/ (Interval_nt a, double b)
    {
      if(b<0){ a = -a; b = -b; }
      else if(b==0) return largest();
      // Now b > 0
      Internal_protector P;
# ifdef __GNUC__
      // Paradoxically, constants should be safe, and this lets the compiler optimize x/2 to x*.5
      if (!__builtin_constant_p(b))
# endif
        b = IA_opacify(b);
      __m128d bb = _mm_set1_pd(b);
      __m128d aa = IA_opacify128(a.simd());
      __m128d r = _mm_div_pd(aa, bb);
      return Interval_nt(IA_opacify128(r));
    }
#endif
};

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
template <bool Protected>
const typename Interval_nt<Protected>::Test_runtime_rounding_modes
Interval_nt<Protected>::tester;
#endif



// Non-documented
// Returns true if the interval is a unique representable double.
template <bool Protected>
inline
bool
fit_in_double (const Interval_nt<Protected> & d, double &r)
{
  bool b = d.is_point();
  if (b)
    r = d.inf();
  return b;
}

// Non-documented
template <bool Protected>
inline
bool
is_singleton (const Interval_nt<Protected> & d)
{
  return d.is_point();
}

// Non-documented
template <bool Protected>
inline
double
magnitude (const Interval_nt<Protected> & d)
{
#ifdef CGAL_USE_SSE2
  const __m128d m = _mm_castsi128_pd (_mm_set1_epi64x (0x7fffffffffffffff));
  __m128d x = _mm_and_pd (d.simd(), m); // { abs(inf), abs(sup) }
  __m128d y = _mm_unpackhi_pd (x, x);
  return _mm_cvtsd_f64 (_mm_max_sd (x, y));
#else
  return (std::max)(CGAL::abs(d.inf()), CGAL::abs(d.sup()));
#endif
}

// Non-documented
template <bool Protected>
inline
double
width (const Interval_nt<Protected> & d)
{
  return d.sup() - d.inf();
}

// Non-documented
template <bool Protected>
inline
double
radius (const Interval_nt<Protected> & d)
{
  return width(d)/2; // This could be improved to avoid overflow.
}

// Non-documented
// This is the relative precision of to_double() (the center of the interval),
// hence we use radius() instead of width().
template <bool Protected>
inline
bool
has_smaller_relative_precision(const Interval_nt<Protected> & d, double prec)
{
  return magnitude(d) == 0 || radius(d) < prec * magnitude(d);
}

// Non-documented
template <bool Protected>
double
relative_precision(const Interval_nt<Protected> & d)
{
  if (magnitude(d) == 0.0)
    return 0.0;
  return radius(d) / magnitude(d);
}


template< bool Protected >
class Is_valid< Interval_nt<Protected> >
  : public CGAL::cpp98::unary_function< Interval_nt<Protected>, bool > {
  public :
    bool operator()( const Interval_nt<Protected>& x ) const {
      return is_valid(-x.inf()) &&
             is_valid(x.sup()) &&
             x.inf() <= x.sup();
    }
};


typedef Interval_nt<false> Interval_nt_advanced;  // for backward-compatibility


// TODO: What about these two guys? Where do they belong to?
template <bool Protected>
struct Min <Interval_nt<Protected> >
    : public CGAL::cpp98::binary_function<Interval_nt<Protected>,
                             Interval_nt<Protected>,
                             Interval_nt<Protected> >
{
    Interval_nt<Protected> operator()( const Interval_nt<Protected>& d,
                                       const Interval_nt<Protected>& e) const
    {
#ifdef CGAL_USE_SSE2
        __m128d x = _mm_min_pd (d.simd(), e.simd());
        // Use _mm_max_sd instead?
        __m128d y = _mm_max_pd (d.simd(), e.simd());
        return Interval_nt<Protected> (_mm_move_sd (x, y));
#else
        return Interval_nt<Protected>(
                -(std::max)(-d.inf(), -e.inf()),
                 (std::min)( d.sup(),  e.sup()));
#endif
    }
};

template <bool Protected>
struct Max <Interval_nt<Protected> >
    : public CGAL::cpp98::binary_function<Interval_nt<Protected>,
                             Interval_nt<Protected>,
                             Interval_nt<Protected> >
{
    Interval_nt<Protected> operator()( const Interval_nt<Protected>& d,
                                       const Interval_nt<Protected>& e) const
    {
#ifdef CGAL_USE_SSE2
        // Use _mm_min_sd instead?
        __m128d x = _mm_min_pd (d.simd(), e.simd());
        __m128d y = _mm_max_pd (d.simd(), e.simd());
        return Interval_nt<Protected> (_mm_move_sd (y, x));
#else
        return Interval_nt<Protected>(
                -(std::min)(-d.inf(), -e.inf()),
                 (std::max)( d.sup(),  e.sup()));
#endif
    }
};

template<bool Protected> inline
Interval_nt<Protected> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Interval_nt<Protected> & x,
const Interval_nt<Protected> & y){
  return CGAL::Min<Interval_nt<Protected> > ()(x,y);
}
template<bool Protected> inline
Interval_nt<Protected> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Interval_nt<Protected> & x,
const Interval_nt<Protected> & y){
  return CGAL::Max<Interval_nt<Protected> > ()(x,y);
}



// TODO : document, when we are OK with the interface.
// - should it allow other number types for the exponent ?
template < bool b >
Interval_nt<b>
ldexp(const Interval_nt<b> &i, int e)
{
  double scale = std::ldexp(1.0, e);
  Interval_nt<b> scale_interval (
                      CGAL_NTS is_finite(scale) ? scale : CGAL_IA_MAX_DOUBLE,
                      scale == 0 ? CGAL_IA_MIN_DOUBLE : scale);
  return i * scale_interval;
}


// We also specialize some corresponding functors returning Uncertain<>.

// TODO: To which concept do these functors belong? Can we remove them??
template < bool b >
struct Equal_to < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x == y; }
};

template < bool b >
struct Not_equal_to < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x != y; }
};

template < bool b >
struct Greater < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x > y; }
};

template < bool b >
struct Less < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x < y; }
};

template < bool b >
struct Greater_equal < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x >= y; }
};

template < bool b >
struct Less_equal < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x <= y; }
};


// As in MP_float.h, the namespace INTERN_INTERVAL_NT contains (now) global
// functions like square or sqrt which would have collided with the new
// global functions from AST/RET
//
// TODO: IMHO, a better solution would be to put the INTERN_MP_FLOAT-functions
//       into the MP_Float-class... But there is surely a reason why this is not
//       the case..?


namespace INTERN_INTERVAL_NT {

  template <bool Protected>
  inline
  double
  to_double (const Interval_nt<Protected> & d)
  {
    return (d.sup() + d.inf()) * 0.5;
    // This may overflow...
  }

  template <bool Protected>
  inline
  std::pair<double, double>
  to_interval (const Interval_nt<Protected> & d)
  {
    return d.pair();
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  sqrt (const Interval_nt<Protected> & d)
  {
    typename Interval_nt<Protected>::Internal_protector P;  // not optimal here.
    // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
    // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
    // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
#ifdef __AVX512F__
    double i = 0;
    if(d.inf() > 0){
      __m128d x = d.simd();
      __m128d m = _mm_set_sd(-0.);
      __m128d y = _mm_xor_pd(x, m);
      // We don't opacify because hopefully a rounded operation is explicit
      // enough that compilers won't mess with it, and it does not care about
      // fesetround.
      __m128d vr = _mm_sqrt_round_sd(y, y, _MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC);
      i = _mm_cvtsd_f64(vr);
      // We could compute the sqrt of d.sup() using _mm_sqrt_pd (same speed as
      // _sd except on broadwell) so it is already in the high part and we can
      // call _mm_sqrt_round_sd(y, x, ...) to merge them directly, but I doubt
      // it helps significantly, it might even hurt by introducing a
      // dependency.
    }
#else
    // TODO: Alternative for computing CGAL_IA_SQRT_DOWN(d.inf()) exactly
    // without changing the rounding mode:
    // - compute x = CGAL_IA_SQRT(d.inf())
    // - compute y = CGAL_IA_SQUARE(x)
    // - if y==d.inf() use x, else use -CGAL_IA_SUB(CGAL_IA_MIN_DOUBLE,x)
    FPU_set_cw(CGAL_FE_DOWNWARD);
    double i = (d.inf() > 0.0) ? CGAL_IA_SQRT(d.inf()) : 0.0;
    FPU_set_cw(CGAL_FE_UPWARD);
#endif
    return Interval_nt<Protected>(i, CGAL_IA_SQRT(d.sup()));
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  square (const Interval_nt<Protected> & d)
  {
    //TODO: SSE version, possibly using abs
    typename Interval_nt<Protected>::Internal_protector P;
    if (d.inf()>=0.0)
        return Interval_nt<Protected>(-CGAL_IA_MUL(-d.inf(), d.inf()),
                                 CGAL_IA_SQUARE(d.sup()));
    if (d.sup()<=0.0)
        return Interval_nt<Protected>(-CGAL_IA_MUL(d.sup(), -d.sup()),
                               CGAL_IA_SQUARE(-d.inf()));
    return Interval_nt<Protected>(0.0, CGAL_IA_SQUARE((std::max)(-d.inf(),
                     d.sup())));
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  abs (const Interval_nt<Protected> & d)
  {
#ifdef CGAL_USE_SSE2
    __m128d a = d.simd();
    __m128d b = (-d).simd();
    __m128d x = _mm_min_pd (a, b);
    __m128d y = _mm_max_pd (a, b);
    __m128d t = _mm_move_sd (y, x);
    __m128d z = _mm_set1_pd(-0.); // +0. would be valid, but I'd rather end up with interval [+0, sup]
    __m128d r = _mm_min_sd(t, z);
    return Interval_nt<Protected> (r);
#else
    if (d.inf() >= 0.0) return d;
    if (d.sup() <= 0.0) return -d;
    return Interval_nt<Protected>(0.0, (std::max)(-d.inf(), d.sup()));
#endif
  }

  template <bool Protected>
  inline
  Uncertain<Sign>
  sign (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 0.0) return POSITIVE;
    if (d.sup() < 0.0) return NEGATIVE;
    if (d.inf() == d.sup()) return ZERO;
    return Uncertain<Sign>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<Comparison_result>
  compare (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
  {
    if (d.inf() > e.sup()) return LARGER;
    if (e.inf() > d.sup()) return SMALLER;
    if (e.inf() == d.sup() && d.inf() == e.sup()) return EQUAL;
    return Uncertain<Comparison_result>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_zero (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 0.0) return false;
    if (d.sup() < 0.0) return false;
    if (d.inf() == d.sup()) return true;
    return Uncertain<bool>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_positive (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 0.0) return true;
    if (d.sup() <= 0.0) return false;
    return Uncertain<bool>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_negative (const Interval_nt<Protected> & d)
  {
    if (d.inf() >= 0.0) return false;
    if (d.sup() < 0.0) return true;
    return Uncertain<bool>::indeterminate();
  }

} // namespace INTERN_INTERVAL_NT


template< bool B > class Real_embeddable_traits< Interval_nt<B> >
  : public INTERN_RET::Real_embeddable_traits_base< Interval_nt<B> , CGAL::Tag_true> {
  public:
    typedef Interval_nt<B>  Type;
  typedef Uncertain<CGAL::Sign> Sign;
  typedef Uncertain<bool> Boolean;
  typedef Uncertain<CGAL::Comparison_result> Comparison_result;

    class Abs
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::abs( x );
        }
    };

    class Sgn
        : public CGAL::cpp98::unary_function< Type, Uncertain< ::CGAL::Sign > > {
      public:
        Uncertain< ::CGAL::Sign > operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::sign( x );
        }
    };

    class Is_positive
      : public CGAL::cpp98::unary_function< Type, Uncertain<bool> > {
      public:
        Uncertain<bool> operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_positive( x );
        }
    };

    class Is_negative
      : public CGAL::cpp98::unary_function< Type, Uncertain<bool> > {
      public:
        Uncertain<bool> operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_negative( x );
        }
    };

    class Compare
      : public CGAL::cpp98::binary_function< Type, Type, Comparison_result > {
      public:
      Comparison_result operator()( const Type& x, const Type& y ) const {
        return INTERN_INTERVAL_NT::compare( x, y );
      }
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
          Comparison_result )
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::to_double( x );
        }
    };

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::to_interval( x );
        }
    };

    class Is_finite
      : public CGAL::cpp98::unary_function< Type, Boolean > {
      public :
        Boolean operator()( const Type& x ) const {
          return CGAL_NTS is_finite( x.inf() ) && CGAL_NTS is_finite( x.sup() );
        }
    };

};

// Algebraic structure traits
template< bool B >
class Algebraic_structure_traits< Interval_nt<B> >
  : public Algebraic_structure_traits_base< Interval_nt<B>,
                                            Field_with_sqrt_tag >  {
  public:
    typedef Interval_nt<B>      Type;
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;
    typedef Uncertain<bool>     Boolean;

    class Is_zero
      : public CGAL::cpp98::unary_function< Type, Boolean > {
      public:
        Boolean operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_zero( x );
        }
    };

    // Specialized just to specify the result type
    class Is_one
      : public CGAL::cpp98::unary_function< Type, Boolean > {
      public:
        Boolean operator()( const Type& x ) const {
          return x == 1;
        }
    };

    class Square
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::square( x );
        }
    };

    class Sqrt
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::sqrt( x );
        }
    };

    struct Is_square
        :public CGAL::cpp98::binary_function<Interval_nt<B>,Interval_nt<B>&,Boolean >
    {
        Boolean operator()(const Interval_nt<B>& x) const {
            return INTERN_INTERVAL_NT::is_positive( x );
        }

        Boolean operator()(
                const Interval_nt<B>& x,
                Interval_nt<B>      & result) const {
            Boolean is_positive = INTERN_INTERVAL_NT::is_positive( x );
            if ( is_positive.inf() == true ){
                typename Algebraic_structure_traits<Interval_nt<B> >::Sqrt sqrt;
                result = sqrt(x);
            }else{
                typename Real_embeddable_traits<Interval_nt<B> >::Abs  abs;
                typename Algebraic_structure_traits<Interval_nt<B> >::Sqrt sqrt;
                result = sqrt(abs(x));
            }
            return is_positive;
        }
    };

  class Divides
    : public CGAL::cpp98::binary_function< Type, Type, Boolean > {
  public:
    Boolean operator()( const Type& x, const Type&) const {
      return ! Is_zero()(x);
    }
    // second operator computing q
    Boolean operator()( const Type& x, const Type& y, Type& q) const {
      if (! Is_zero()(x) )
        q  = y/x ;
      return Boolean(true);
    }
  };
};


// COERCION_TRAITS BEGIN
template < class A, class B , int > struct Coercion_traits_for_level;
template < class A, class B, class C> struct Coercion_traits_interval_nt;

template<class A ,bool P >
struct Coercion_traits_for_level<A,Interval_nt<P>,CTL_INTERVAL>
    :public Coercion_traits_interval_nt<A,Interval_nt<P>,
            typename Real_embeddable_traits<A>::Is_real_embeddable>{};

template<class A , bool P>
struct Coercion_traits_for_level<Interval_nt<P>,A,CTL_INTERVAL>
    :public Coercion_traits_for_level<A,Interval_nt<P>, CTL_INTERVAL>{};

template<class A , bool P >
struct Coercion_traits_interval_nt<A, Interval_nt<P>,Tag_false>
    :public Coercion_traits_for_level<A,Interval_nt<P>,0>{};

template<class A , bool P>
struct Coercion_traits_interval_nt<A, Interval_nt<P>, Tag_true>{
    typedef Tag_true Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef Interval_nt<P> Type;
    struct Cast {
        typedef Interval_nt<P> result_type;
        Interval_nt<P> inline operator()(const Interval_nt<P>& x ) const {
            return x;
        }
        Interval_nt<P> inline operator()(const A& x ) const {
            return typename Real_embeddable_traits<A>::To_interval()(x);
        }
    };
};

// COERCION_TRAITS END

template< bool B >
class Interval_traits< Interval_nt<B> >
  : public internal::Interval_traits_base< Interval_nt<B> >  {
public:
  typedef Interval_traits<Interval_nt<B> > Self;
  typedef Interval_nt<B> Interval;
  typedef double Bound;
  typedef CGAL::Tag_false With_empty_interval;
  typedef CGAL::Tag_true  Is_interval;

 struct Construct :public CGAL::cpp98::binary_function<Bound,Bound,Interval>{
    Interval operator()( const Bound& l,const Bound& r) const {
      CGAL_precondition( l < r );
      return Interval(l,r);
    }
  };

  struct Lower :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.inf();
    }
  };

  struct Upper :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.sup();
    }
  };

  struct Width :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return width(a);
    }
  };

  struct Median :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return (Lower()(a)+Upper()(a))/2.0;
    }
  };

  struct Norm :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return magnitude(a);
    }
  };

  struct Singleton :public CGAL::cpp98::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return a.is_point();
    }
  };

  struct Zero_in :public CGAL::cpp98::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return Lower()(a) <= 0  &&  0 <= Upper()(a);
    }
  };

  struct In :public CGAL::cpp98::binary_function<Bound,Interval,bool>{
    bool operator()( Bound x, const Interval& a ) const {
      return Lower()(a) <= x && x <= Upper()(a);
    }
  };

  struct Equal :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.is_same(b);
    }
  };

  struct Overlap :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.do_overlap(b);
    }
  };

  struct Subset :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Lower()(b) <= Lower()(a) && Upper()(a) <= Upper()(b) ;
    }
  };

  struct Proper_subset :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Subset()(a,b) && ! Equal()(a,b);
    }
  };

  struct Hull :public CGAL::cpp98::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
#ifdef CGAL_USE_SSE2
      return Interval(_mm_max_pd(a.simd(), b.simd()));
#else
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      return Interval(
             -max BOOST_PREVENT_MACRO_SUBSTITUTION (-a.inf(),-b.inf()),
              max BOOST_PREVENT_MACRO_SUBSTITUTION ( a.sup(), b.sup()));
#endif
    }
  };


//  struct Empty is Null_functor

  struct Intersection :public CGAL::cpp98::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      Bound l(max BOOST_PREVENT_MACRO_SUBSTITUTION (Lower()(a),Lower()(b)));
      Bound u(min BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(a),Upper()(b)));
      if(u < l ) throw Exception_intersection_is_empty();
      return Construct()(l,u);
    }
  };
};

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<bool b> struct NumTraits<CGAL::Interval_nt<b> >
  {
    typedef CGAL::Interval_nt<b> Real;
    typedef CGAL::Interval_nt<b> NonInteger;
    typedef CGAL::Interval_nt<b> Nested;
    typedef double Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }
    static inline Real highest() { return Real((std::numeric_limits<double>::max)(), std::numeric_limits<double>::infinity()); }
    static inline Real lowest() { return Real(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::lowest()); }

    // Costs could depend on b.
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 0,
      ReadCost = 2,
      AddCost = 2,
      MulCost = 10
    };
  };

  namespace internal {
    template<class> struct significant_decimals_impl;
    template<bool b>
      struct significant_decimals_impl<CGAL::Interval_nt<b> >
      : significant_decimals_impl<typename CGAL::Interval_nt<b>::value_type> { };

    // Without this, when computing some decompositions for a matrix of
    // intervals, Eigen looks for the largest element in a column (for
    // instance). There may easily be 2 equal, slightly imprecise numbers that
    // could equally well be used as pivots, but Eigen ends up spuriously
    // throwing in the comparison between them. So we provide a different
    // strategy for picking the pivot.
    template<typename> struct scalar_score_coeff_op;
    template<bool b> struct scalar_score_coeff_op<CGAL::Interval_nt<b> > {
      // If all coeffs can be 0, it is essential to designate as the best one
      // that can be non-zero and has a non-zero score, if there is one.
      struct result_type : boost::totally_ordered1<result_type> {
        CGAL::Interval_nt<b> i;
        result_type():i(){}
        result_type(CGAL::Interval_nt<b> j):i(j){}
        friend bool operator<(result_type x, result_type y){
          if(x.i.inf()==0){
            if(y.i.inf()==0)return x.i.sup()<y.i.sup(); // [0,0]<[0,1]
            else return true; // [0,*]<[1,*]
          }
#if 0
          // The following is already handled by the general formula below
          if(y.i.inf()==0)return false; // [0,*]<[1,*]
#endif
          // Both numbers are guaranteed non-zero. With double people usually
          // pick the biggest number. Here we choose the tightest interval.
          // This is purely heuristic, it doesn't matter much if overflow makes
          // us do random choices.
          // Best is largest inf/sup (ideally 1)
          // Risk of {over,under}flow
          return x.i.inf()*y.i.sup() < y.i.inf()*x.i.sup();
        }
        // Only used as: if(max==Score(0))
        friend bool operator==(result_type x, result_type y){
          // Throw if we don't know if the max coeff is 0
          return x.i == y.i;
        }
      };
      result_type operator()(CGAL::Interval_nt<b> const&x)const{return abs(x);}
    };
    template<typename> struct functor_traits;
    template<bool b> struct functor_traits<scalar_score_coeff_op<CGAL::Interval_nt<b> > >
    {
      enum {
        Cost = 10,
        PacketAccess = false
      };
    };
  }
}

#undef CGAL_CST_TRUE

#endif // CGAL_INTERVAL_NT_H
