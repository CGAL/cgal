// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_FILTERED_EXTENDED_HOMOGENEOUS_H
#define CGAL_FILTERED_EXTENDED_HOMOGENEOUS_H

#include <CGAL/license/Nef_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/number_utils.h>
#include <CGAL/tss.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 5
#include <CGAL/Nef_2/debug.h>

#define REDUCE_INTERSECTION_POINTS
//#define KERNEL_ANALYSIS
//#define KERNEL_CHECK

#ifdef  KERNEL_CHECK
#include <CGAL/Extended_homogeneous.h>
#define CHECK(c1,c2) CGAL_assertion((c1) == (c2));
#define PRINT_CHECK_ENABLED std::cout << "kernel check enabled!\n"
#else
#define CHECK(c1,c2)
#define PRINT_CHECK_ENABLED std::cout << "no kernel check!\n"
#endif

#ifdef KERNEL_ANALYSIS
#define DEFCOUNTER(c) \
  static int c##_total=0; static int c##_exception=0;
#define INCTOTAL(c) c##_total++
#define INCEXCEPTION(c) c##_exception++
#define PRINT_STATISTICS(c) \
std::cout << #c##" " << c##_exception << "/" << c##_total << std::endl
#else
#define DEFCOUNTER(c)
#define INCTOTAL(c)
#define INCEXCEPTION(c)
#define PRINT_STATISTICS(c)
#endif

namespace CGAL {

template <typename RT>
class SPolynomial {
  RT _m,_n;
public:
  SPolynomial() : _m(),_n() {}
  SPolynomial(const RT& m, const RT& n) : _m(m),_n(n) {}
  SPolynomial(const RT& n) : _m(),_n(n) {}
  SPolynomial(const SPolynomial<RT>& p) : _m(p._m),_n(p._n) {}
  SPolynomial<RT>& operator=(const SPolynomial<RT>& p)
  { _m=p._m; _n=p._n; return *this; }

  const RT& m() const { return _m; }
  const RT& n() const { return _n; }
  void negate() { _m=-_m; _n=-_n; }

  SPolynomial<RT> operator*(const RT& c) const
  { return SPolynomial<RT>(c*_m,c*_n); }
  SPolynomial<RT> operator+(const SPolynomial<RT>& p) const
  { return SPolynomial<RT>(_m+p._m,_n+p._n); }
  SPolynomial<RT> operator-(const SPolynomial<RT>& p) const
  { return SPolynomial<RT>(_m-p._m,_n-p._n); }
  SPolynomial<RT> operator-() const
  { return SPolynomial<RT>(-_m,-_n); }
  void operator /= (const RT& c)
  { _m /= c; _n /= c; }
  const RT& operator[](int i) { return (i%2 ? _n : _m); }
  const RT& operator[](int i) const { return (i%2 ? _n : _m); }
  bool is_zero() const { return (_m==0 && _n==0); }
  int sign() const
  { if ( _m != 0 ) return CGAL_NTS sign(_m);
    return CGAL_NTS sign(_n);
  }


  // only for visualization:

  static RT& R()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(RT,R_,0);
  }

  static void set_R(const RT& r) { R() = r; }
  RT eval_at(const RT& r) const { return _m*r+_n; }
  RT eval_at_R() const { return _m*R()+_n; }

};


template <typename RT>
int sign(const SPolynomial<RT>& p)
{
  return p.sign();
}


template <typename RT>
bool operator==(const SPolynomial<RT>& p1, const SPolynomial<RT>& p2)
{ return (p1-p2).is_zero(); }

template <typename RT>
bool operator>(const SPolynomial<RT>& p1, const SPolynomial<RT>& p2)
{ return (p1-p2).sign()>0; }

template <typename RT>
bool operator<(const SPolynomial<RT>& p1, const SPolynomial<RT>& p2)
{ return (p1-p2).sign()<0; }

template <typename RT>
bool operator<(int i, const SPolynomial<RT>& p2)
{
  SPolynomial<RT> p1(i);
  return (p1-p2).sign()<0;
}


template <typename RT>
bool operator<(const SPolynomial<RT>& p1, int i)
{
  SPolynomial<RT> p2(i);
  return (p1-p2).sign()<0;
}

template <class RT>
inline double to_double(const SPolynomial<RT>& p)
{ return (CGAL::to_double(p.eval_at(SPolynomial<RT>::R()))); }

template <class RT>
std::ostream& operator<<(std::ostream& os, const SPolynomial<RT>& p)
{
  switch( get_mode(os) ) {
    case CGAL::IO::ASCII :
      os << p.m() << " " << p.n(); break;
    case CGAL::IO::BINARY :
      CGAL::write(os,p.m());CGAL::write(os,p.n()); break;
    default:
      if ( p.m() == 0 ) os<<"["<<p.n()<<"]";
      else os<<"["<<p.m()<<" R + "<<p.n()<<"]";
  }
  return os;
}
template <class RT>
std::istream& operator>>(std::istream& is, SPolynomial<RT>& p)
{ RT m,n;
  switch( get_mode(is) ){
    case CGAL::IO::ASCII :
      is >> m >> n; p = SPolynomial<RT>(m,n); break;
    case CGAL::IO::BINARY :
      CGAL::read(is,m);CGAL::read(is,n);break;
    default:
    CGAL_error_msg("\nStream must be in ascii or binary mode\n");
      break;
  }
  return is;
}

template <class RT> /*CGAL_KERNEL_INLINE*/
CGAL::io_Operator io_tag(const SPolynomial<RT>&)
{ return CGAL::io_Operator(); }


template <typename RT>
class SQuotient {
  SPolynomial<RT> _p;
  RT              _n;
public:
  SQuotient() : _p(),_n() {}
  SQuotient(const SPolynomial<RT>& p, const RT& n) : _p(p),_n(n) {}
  SQuotient(const SQuotient<RT>& p) : _p(p._p),_n(p._n) {}
  SQuotient<RT>& operator=(const SQuotient<RT>& p)
  { _p=p._p; _n=p._n; return *this; }
  const SPolynomial<RT>& numerator() const { return _p; }
  const RT&              denominator() const { return _n; }
};

template <class RT>
inline double to_double(const SQuotient<RT>& q)
{ return (CGAL::to_double(q.numerator().eval_at_R())/
          CGAL::to_double(q.denominator())); }


template <typename RT> class Extended_point;
template <typename RT> class Extended_point_rep;

template <typename RT>
class Extended_point_rep {
  friend class Extended_point<RT>;
  SPolynomial<RT> x_,y_; RT w_;
  typedef Interval_nt_advanced DT;
  DT mxd,myd,nxd,nyd,wd;
public:
  Extended_point_rep(const RT& x, const RT& y, const RT& w) :
    x_(x),y_(y),w_(w)
  { CGAL_assertion_msg(w!=0,"denominator is zero.");
    nxd=CGAL::to_interval(x);
    nyd=CGAL::to_interval(y);
    wd=CGAL::to_interval(w);
    mxd=myd=0;
  }

  Extended_point_rep(const SPolynomial<RT>& x,
                     const SPolynomial<RT>& y,
                     const RT& w) : x_(x),y_(y),w_(w)
  { CGAL_assertion_msg(w!=0,"denominator is zero.");
    mxd=CGAL::to_interval(x.m());
    myd=CGAL::to_interval(y.m());
    nxd=CGAL::to_interval(x.n());
    nyd=CGAL::to_interval(y.n());
    wd=CGAL::to_interval(w);
  }

  Extended_point_rep() : x_(),y_(),w_() {}
  ~Extended_point_rep() {}
  void negate()
  { x_ = -x_; y_ = -y_; w_ = -w_;
    mxd = -mxd; myd = -myd; nxd = -nxd; nyd = -nyd; wd = -wd; }

};

template <typename RT_>
class Extended_point : public Handle_for< Extended_point_rep<RT_> > {
  typedef Extended_point_rep<RT_> Rep;
  typedef Handle_for< Rep >       Base;

  using Base::ptr;

public:
  typedef typename Rep::DT DT;
  typedef RT_ RT;
  typedef SPolynomial<RT>  SP;

  Extended_point() : Base( Rep() ) {}

  Extended_point(const RT& x, const RT& y, const RT& w) :
    Base( Rep(x,y,w) )
  { if (w < 0) this->ptr()->negate(); }

  Extended_point(const SPolynomial<RT>& x,
                 const SPolynomial<RT>& y,
                 const RT& w) : Base( Rep(x,y,w) )
  { if (w < 0) this->ptr()->negate(); }

  Extended_point(const RT& mx, const RT& nx,
                 const RT& my, const RT& ny, const RT& w) :
    Base( Rep(SP(mx,nx), SP(my,ny), w) )
  { if (w < 0) this->ptr()->negate(); }

  Extended_point(const Extended_point<RT>& p) : Base(p) {}
  ~Extended_point() {}

  Extended_point& operator=(const Extended_point<RT>& p)
  { Base::operator=(p); return *this; }

  const RT& mx() const { return this->ptr()->x_.m(); }
  const RT& nx() const { return this->ptr()->x_.n(); }
  const RT& my() const { return this->ptr()->y_.m(); }
  const RT& ny() const { return this->ptr()->y_.n(); }
  const RT& hw()  const { return this->ptr()->w_; }
  const DT& mxD() const { return this->ptr()->mxd; }
  const DT& nxD() const { return this->ptr()->nxd; }
  const DT& myD() const { return this->ptr()->myd; }
  const DT& nyD() const { return this->ptr()->nyd; }
  const DT& hwD() const { return this->ptr()->wd; }

  SQuotient<RT> x() const
  { return SQuotient<RT>(this->ptr()->x_, this->ptr()->w_); }
  SQuotient<RT> y() const
  { return SQuotient<RT>(this->ptr()->y_, this->ptr()->w_); }

  const SPolynomial<RT> hx() const { return this->ptr()->x_; }
  const SPolynomial<RT> hy() const { return this->ptr()->y_; }

  bool is_standard() const { return (mx()==0)&&(my()==0); }
  Extended_point<RT> opposite() const
  { return Extended_point<RT>(-mx(),nx(),-my(),ny(),this->w()); }

#ifdef KERNEL_CHECK
typedef CGAL::Extended_homogeneous<RT_> CheckKernel;
typedef typename CheckKernel::Point_2   CheckPoint;
typedef typename CheckKernel::RT        CheckRT;

CheckRT convert(const CGAL::SPolynomial<RT_>& p) const
{ return CheckRT(p.n(),p.m()); }
CheckRT convert(const RT_& t) const
{ return CheckRT(t); }
CheckPoint checkrep() const
{ return CheckPoint(convert(hx()),convert(hy()),convert(w())); }

#endif // KERNEL_CHECK


};

template <class RT>
std::ostream& operator<<(std::ostream& os, const Extended_point<RT>& p)
{ switch( get_mode(os) ) {
    case CGAL::IO::ASCII :
      os << p.hx() << " " << p.hy() << " " << p.hw(); break;
    case CGAL::IO::BINARY :
      CGAL::write(os,p.hx());CGAL::write(os,p.hy());
      CGAL::write(os,p.hw()); break;
    default:
      os << "(" << p.hx() << "," << p.hy() << "," << p.hw() << ")";
#if 0
      os << "((" << CGAL::to_double(p.nx())/CGAL::to_double(p.hw()) << ","
         << CGAL::to_double(p.ny())/CGAL::to_double(p.hw()) << "))";
#endif
  }
  return os;
}
template <class RT>
std::istream& operator>>(std::istream& is, Extended_point<RT>& p)
{ SPolynomial<RT> x,y; RT w;
  switch( get_mode(is) ){
    case CGAL::IO::ASCII :
      is >> x >> y >> w; break;
    case CGAL::IO::BINARY :
      CGAL::read(is,x);CGAL::read(is,y);CGAL::read(is,w); break;
    default:
    CGAL_error_msg("\nStream must be in ascii or binary mode\n");
      break;
  }
  p = Extended_point<RT>(x,y,w);
  return is;
}


template <typename NT> inline
int orientation_coeff2(const NT& mx1, const NT& /*nx1*/,
                       const NT& my1, const NT& /*ny1*/, const NT& w1,
                       const NT& mx2, const NT& /*nx2*/,
                       const NT& my2, const NT& /*ny2*/, const NT& w2,
                       const NT& mx3, const NT& /*nx3*/,
                       const NT& my3, const NT& /*ny3*/, const NT& w3)
{
  NT coeff2 = mx1*w3*my2-mx1*w2*my3+mx3*w2*my1-
              mx2*w3*my1-mx3*w1*my2+mx2*w1*my3;
  return CGAL_NTS sign(coeff2);
}

template <typename NT> inline
int orientation_coeff1(const NT& mx1, const NT& nx1,
                       const NT& my1, const NT& ny1, const NT& w1,
                       const NT& mx2, const NT& nx2,
                       const NT& my2, const NT& ny2, const NT& w2,
                       const NT& mx3, const NT& nx3,
                       const NT& my3, const NT& ny3, const NT& w3)
{
  NT coeff1 = mx1*w3*ny2-mx1*w2*ny3+nx1*w3*my2-mx2*w3*ny1-
              nx1*w2*my3+mx2*w1*ny3-nx2*w3*my1+mx3*w2*ny1+
              nx2*w1*my3-mx3*w1*ny2+nx3*w2*my1-nx3*w1*my2;
  return CGAL_NTS sign(coeff1);
}

template <typename NT> inline
int orientation_coeff0(const NT& /*mx1*/, const NT& nx1,
                       const NT& /*my1*/, const NT& ny1, const NT& w1,
                       const NT& /*mx2*/, const NT& nx2,
                       const NT& /*my2*/, const NT& ny2, const NT& w2,
                       const NT& /*mx3*/, const NT& nx3,
                       const NT& /*my3*/, const NT& ny3, const NT& w3)
{
  NT coeff0 = -nx2*w3*ny1+nx1*w3*ny2+nx2*w1*ny3-
               nx1*w2*ny3+nx3*w2*ny1-nx3*w1*ny2;
  return CGAL_NTS sign(coeff0);
}

DEFCOUNTER(or0)
DEFCOUNTER(or1)
DEFCOUNTER(or2)
template <typename RT>
int orientation(const Extended_point<RT>& p1,
                const Extended_point<RT>& p2,
                const Extended_point<RT>& p3)
{ CGAL_NEF_TRACEN("orientation "<<p1<<p2<<p3);
  int res;
  try { INCTOTAL(or2); Protect_FPU_rounding<true> Protection;
    res = orientation_coeff2(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                             p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                             p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(or2);
    res = orientation_coeff2(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                             p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                             p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(or1); Protect_FPU_rounding<true> Protection;
    res = orientation_coeff1(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                             p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                             p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(or1);
    res = orientation_coeff1(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                             p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                             p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(or0); Protect_FPU_rounding<true> Protection;
    res = orientation_coeff0(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                             p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                             p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(or0);
    res = orientation_coeff0(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                             p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                             p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw());
  }
  return res;
}

template <typename NT>
inline
int compare_expr(const NT& n1, const NT& d1,
                 const NT& n2, const NT& d2)
{ return CGAL_NTS sign( n1*d2 - n2*d1 ); }

DEFCOUNTER(cmpx0)
DEFCOUNTER(cmpx1)

template <typename RT>
int compare_x(const Extended_point<RT>& p1,
              const Extended_point<RT>& p2)
{
  int res;
  try { INCTOTAL(cmpx1); Protect_FPU_rounding<true> Protection;
    res = compare_expr(p1.mxD(),p1.hwD(),p2.mxD(),p2.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmpx1);
    res = compare_expr(p1.mx(),p1.hw(),p2.mx(),p2.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(cmpx0); Protect_FPU_rounding<true> Protection;
    res = compare_expr(p1.nxD(),p1.hwD(),p2.nxD(),p2.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmpx0);
    res = compare_expr(p1.nx(),p1.hw(),p2.nx(),p2.hw());
  }
  return res;
}

DEFCOUNTER(cmpy0)
DEFCOUNTER(cmpy1)

template <typename RT>
int compare_y(const Extended_point<RT>& p1,
              const Extended_point<RT>& p2)
{
  int res;
  try { INCTOTAL(cmpy1); Protect_FPU_rounding<true> Protection;
    res = compare_expr(p1.myD(),p1.hwD(),p2.myD(),p2.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmpy1);
    res = compare_expr(p1.my(),p1.hw(),p2.my(),p2.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(cmpy0); Protect_FPU_rounding<true> Protection;
    res = compare_expr(p1.nyD(),p1.hwD(),p2.nyD(),p2.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmpy0);
    res = compare_expr(p1.ny(),p1.hw(),p2.ny(),p2.hw());
  }
  return res;
}


template <typename RT>
inline
int compare_xy(const Extended_point<RT>& p1,
               const Extended_point<RT>& p2)
{ int c1 = compare_x(p1,p2);
  if ( c1 != 0 ) return c1;
  else return compare_y(p1,p2);
}

template <typename RT>
inline
bool strictly_ordered_along_line(const Extended_point<RT>& p1,
                                 const Extended_point<RT>& p2,
                                 const Extended_point<RT>& p3)
{ return ( orientation(p1,p2,p3) == 0 ) &&
         ( compare_xy(p1,p2) * compare_xy(p2,p3) == 1 );
}

template <typename RT>
inline bool operator==(const Extended_point<RT>& p1,
                       const Extended_point<RT>& p2)
{ CHECK(bool(compare_xy(p1,p2) == 0),p1.checkrep()==p2.checkrep())
  return (p1.identical(p2) || compare_xy(p1,p2) == 0); }

template <typename RT>
inline bool operator!=(const Extended_point<RT>& p1,
                       const Extended_point<RT>& p2)
{ return !(p1==p2); }


template <typename NT>
inline
int cmppd_coeff2(const NT& mx1, const NT& /*nx1*/,
                 const NT& my1, const NT& /*ny1*/, const NT& w1,
                 const NT& mx2, const NT& /*nx2*/,
                 const NT& my2, const NT& /*ny2*/, const NT& w2,
                 const NT& mx3, const NT& /*nx3*/,
                 const NT& my3, const NT& /*ny3*/, const NT& w3,
                 const NT& mx4, const NT& /*nx4*/,
                 const NT& my4, const NT& /*ny4*/, const NT& w4)
{
  NT w1Q(w1*w1), w2Q(w2*w2), w3Q(w3*w3), w4Q(w4*w4);
  NT w1w2Q(w1Q*w2Q), w3w4Q(w3Q*w4Q), two(2);
  NT coeff2 =    w3w4Q * w2Q *mx1*mx1-
                 two* w3w4Q  *w2*mx1*w1*mx2+
                 w3w4Q * w1Q *mx2*mx2+
                 w3w4Q * w2Q *my1*my1-
                 two* w3w4Q  *w2*my1*w1*my2+
                 w3w4Q * w1Q *my2*my2-
                 w1w2Q * w4Q *mx3*mx3+
                 two* w1w2Q  *w4*mx3*w3*mx4-
                 w1w2Q * w3Q *mx4*mx4-
                 w1w2Q * w4Q *my3*my3+
                 two* w1w2Q  *w4*my3*w3*my4-
                 w1w2Q * w3Q *my4*my4;
  return CGAL_NTS sign(coeff2);
}


template <typename NT>
inline
int cmppd_coeff1(const NT& mx1, const NT& nx1,
                 const NT& my1, const NT& ny1, const NT& w1,
                 const NT& mx2, const NT& nx2,
                 const NT& my2, const NT& ny2, const NT& w2,
                 const NT& mx3, const NT& nx3,
                 const NT& my3, const NT& ny3, const NT& w3,
                 const NT& mx4, const NT& nx4,
                 const NT& my4, const NT& ny4, const NT& w4)
{
  NT w1Q(w1*w1), w2Q(w2*w2), w3Q(w3*w3), w4Q(w4*w4);
  NT w1w2Q(w1Q*w2Q), w3w4Q(w3Q*w4Q), two(2);
  NT coeff1 = two * (w3w4Q * w1Q * mx2*nx2-
                     w3w4Q * w2*my1*w1*ny2+
                     w3w4Q * w1Q * my2*ny2+
                     w1w2Q * w4*nx3*w3*mx4-
                     w1w2Q * w4Q *mx3*nx3+
                     w3w4Q * w2Q *mx1*nx1-
                     w3w4Q * w2*mx1*w1*nx2-
                     w3w4Q * w2*nx1*w1*mx2-
                     w3w4Q * w2*ny1*w1*my2-
                     w1w2Q * w4Q *my3*ny3+
                     w1w2Q * w4*my3*w3*ny4+
                     w1w2Q * w4*ny3*w3*my4+
                     w3w4Q * w2Q *my1*ny1-
                     w1w2Q * w3Q *my4*ny4+
                     w1w2Q * w4*mx3*w3*nx4-
                     w1w2Q * w3Q *mx4*nx4);
  return CGAL_NTS sign(coeff1);
}

template <typename NT>
inline
int cmppd_coeff0(const NT& /*mx1*/, const NT& nx1,
                 const NT& /*my1*/, const NT& ny1, const NT& w1,
                 const NT& /*mx2*/, const NT& nx2,
                 const NT& /*my2*/, const NT& ny2, const NT& w2,
                 const NT& /*mx3*/, const NT& nx3,
                 const NT& /*my3*/, const NT& ny3, const NT& w3,
                 const NT& /*mx4*/, const NT& nx4,
                 const NT& /*my4*/, const NT& ny4, const NT& w4)
{
  NT w1Q(w1*w1), w2Q(w2*w2), w3Q(w3*w3), w4Q(w4*w4);
  NT w1w2Q(w1Q*w2Q), w3w4Q(w3Q*w4Q), two(2);
  NT coeff0 = w3w4Q * (w1Q * ( nx2*nx2 + ny2*ny2 ) +
                       w2Q * ( ny1*ny1 + nx1*nx1 )) -
              w1w2Q * (w4Q * ( nx3*nx3 + ny3*ny3 ) +
                       w3Q * ( nx4*nx4 + ny4*ny4 )) +
              two* (- w3w4Q * (w2*nx1*w1*nx2 + w2*ny1*w1*ny2)
                    + w1w2Q * (w4*ny3*w3*ny4 + w4*nx3*w3*nx4));
  return CGAL_NTS sign(coeff0);
}

DEFCOUNTER(cmppd2)
DEFCOUNTER(cmppd1)
DEFCOUNTER(cmppd0)

// leghth.mws
template <typename RT>
int compare_pair_dist(
  const Extended_point<RT>& p1, const Extended_point<RT>& p2,
  const Extended_point<RT>& p3, const Extended_point<RT>& p4)
{
  int res;
  try { INCTOTAL(cmppd2); Protect_FPU_rounding<true> Protection;
    res = cmppd_coeff2(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                       p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                       p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                       p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmppd2);
    res = cmppd_coeff2(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                       p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                       p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                       p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(cmppd1); Protect_FPU_rounding<true> Protection;
    res = cmppd_coeff1(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                       p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                       p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                       p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmppd1);
    res = cmppd_coeff1(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                       p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                       p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                       p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(cmppd0); Protect_FPU_rounding<true> Protection;
    res = cmppd_coeff0(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                       p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                       p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                       p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  }
  catch (Uncertain_conversion_exception&) { INCEXCEPTION(cmppd0);
    res = cmppd_coeff0(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                       p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                       p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                       p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  return res;
}


template <typename RT>
class Extended_segment {
  Extended_point<RT> _p1,_p2;
public:
  Extended_segment() : _p1(),_p2() {}
  Extended_segment(const Extended_point<RT>& p1,
                   const Extended_point<RT>& p2) :
    _p1(p1), _p2(p2) {}
  Extended_segment(const Extended_segment<RT>& s) :
    _p1(s._p1), _p2(s._p2) {}
  Extended_segment<RT>& operator=(const Extended_segment<RT>& s)
  { _p1 = s._p1; _p2 = s._p2; return *this; }

  const Extended_point<RT>& source() const { return _p1; }
  const Extended_point<RT>& target() const { return _p2; }

  void line_equation(RT& a, RT& b, SPolynomial<RT>& c) const;
};

template <class RT>
std::ostream& operator<<(std::ostream& os, const Extended_segment<RT>& s)
{ os << s.source() << s.target(); return os; }
template <class RT>
std::istream& operator>>(std::istream& is, Extended_segment<RT>& s)
{ Extended_point<RT> p1,p2;
  is >> p1 >> p2; s=Extended_segment<RT>(p1,p2); return is; }

template <typename RT>
void Extended_segment<RT>::
line_equation(RT& a, RT& b, SPolynomial<RT>& c) const
{
  bool sstandard = _p1.is_standard();
  bool tstandard = _p2.is_standard();
  if (sstandard && tstandard) {
    a = _p1.ny()*_p2.hw() - _p2.ny()*_p1.hw();
    b = _p1.hw()*_p2.nx() - _p2.hw()*_p1.nx();
    c = SPolynomial<RT>(_p1.nx()*_p2.ny() - _p2.nx()*_p1.ny());
    return;

  }
  Extended_point<RT> p;
  bool correct_orientation=true;
  if (!sstandard && !tstandard) {
    bool x_equal = (_p1.hx()*_p2.hw() - _p2.hx()*_p1.hw()).is_zero();
    bool y_equal = (_p1.hy()*_p2.hw() - _p2.hy()*_p1.hw()).is_zero();
    if (x_equal && CGAL_NTS abs(_p1.mx())==_p1.hw() && _p1.nx()==0 )
    { int dy = (_p2.hy()-_p1.hy()).sign();
      a=-dy; b=0; c = SPolynomial<RT>(dy*_p1.hx().sign(),0); return; }
    if (y_equal && CGAL_NTS abs(_p1.my())==_p1.hw() && _p1.ny()==0 )
    { int dx = (_p2.hx()-_p1.hx()).sign();
      a=0; b=dx; c = SPolynomial<RT>(-dx*_p1.hy().sign(),0); return; }
    p = _p2; // evaluation according to mixed case

  }
  else if (sstandard && !tstandard)
  { p = _p2; }
  else if (!sstandard && tstandard)
  { p = _p1; correct_orientation=false; }
  RT w = p.hw();
  RT ci;
  if ( correct_orientation ) {
    a = -p.my();
    b =  p.mx();
    ci = (p.nx()*p.my()-p.ny()*p.mx())/w;
  } else {
    a =  p.my();
    b = -p.mx();
    ci = (p.ny()*p.mx()-p.nx()*p.my())/w;
  }
  c = SPolynomial<RT>(ci);



}

template <typename RT>
Extended_point<RT> intersection(
  const Extended_segment<RT>& s1, const Extended_segment<RT>& s2)
{
  RT a1,b1,a2,b2;
  SPolynomial<RT> c1,c2;
  s1.line_equation(a1,b1,c1);
  s2.line_equation(a2,b2,c2);
  SPolynomial<RT> x = c2*b1 - c1*b2;
  SPolynomial<RT> y = c1*a2 - c2*a1;
  RT w = a1*b2 - a2*b1; CGAL_assertion(w!=0);
  #ifdef REDUCE_INTERSECTION_POINTS
  RT xgcd,ygcd;
  if ( x.m() == RT(0) )  xgcd = ( x.n() == 0 ? RT(1) : x.n() );
  else /* != 0 */    xgcd = ( x.n() == 0 ? x.m() : CGAL_NTS gcd(x.m(),x.n()) );
  if ( y.m() == RT(0) )  ygcd = ( y.n() == 0 ? RT(1) : y.n() );
  else /* != 0 */    ygcd = ( y.n() == 0 ? y.m() : CGAL_NTS gcd(y.m(),y.n()) );
  RT d = CGAL_NTS gcd(w,CGAL_NTS gcd(xgcd,ygcd));
  x /= d;
  y /= d;
  w /= d;
  #endif // REDUCE_INTERSECTION_POINTS
  return Extended_point<RT>(x,y,w);
}

template <typename RT>
inline
int orientation(const Extended_segment<RT>& s, const Extended_point<RT>& p)
{ return orientation(s.source(),s.target(),p); }

template <typename RT>
inline
bool is_degenerate(const Extended_segment<RT>& s)
{ return s.source()==s.target(); }

template <typename RT>
inline
bool contains(const Extended_segment<RT>& s,
              const Extended_point<RT>& p)
{ int p_rel_source = compare_xy(p,s.source());
  int p_rel_target = compare_xy(p,s.target());
  return ( orientation(s,p) == 0 ) &&
    ( ( p_rel_source >= 0 && p_rel_target <= 0 ) ||
      ( p_rel_source <= 0 && p_rel_target >= 0 ) );
}


template <typename RT>
class Extended_direction {
  Extended_point<RT> _p1,_p2;
public:
  Extended_direction() : _p1(),_p2() {}
  Extended_direction(const Extended_direction<RT>& d) :
    _p1(d._p1),_p2(d._p2) {}
  Extended_direction<RT>& operator=(const Extended_direction<RT>& d)
  { _p1 = d._p1; _p2 = d._p2; return *this; }

  Extended_direction(const Extended_point<RT>& p1,
                     const Extended_point<RT>& p2) :
    _p1(p1),_p2(p2) {}

  Extended_direction(const RT& x, const RT& y) :
    _p1(0,0,1),_p2(x,y,1) {}

  const Extended_point<RT>& p1() const { return _p1; }
  const Extended_point<RT>& p2() const { return _p2; }
  int dx_sign() const
  { return (_p2.hx()*_p1.hw()-_p1.hx()*_p2.hw()).sign(); }
  int dy_sign() const
  { return (_p2.hy()*_p1.hw()-_p1.hy()*_p2.hw()).sign(); }
};

template <class RT>
std::ostream& operator<<(std::ostream& os, const Extended_direction<RT>& d)
{ os << d.p1() << "," << d.p2();
  return os; }
template <class RT>
std::istream& operator>>(std::istream& is, Extended_direction<RT>& d)
{ Extended_point<RT> x,y;
  is >> x >> y; d = Extended_direction<RT>(x,y);
  return is; }


template <typename NT>
inline
int coeff2_dor(const NT& mx1, const NT& /*nx1*/,
               const NT& my1, const NT& /*ny1*/, const NT& w1,
               const NT& mx2, const NT& /*nx2*/,
               const NT& my2, const NT& /*ny2*/, const NT& w2,
               const NT& mx3, const NT& /*nx3*/,
               const NT& my3, const NT& /*ny3*/, const NT& w3,
               const NT& mx4, const NT& /*nx4*/,
               const NT& my4, const NT& /*ny4*/, const NT& w4)
{
  NT coeff2 = w1*mx2*w3*my4-w1*mx2*w4*my3-w2*mx1*w3*my4+w2*mx1*w4*my3-
              w1*my2*w3*mx4+w1*my2*w4*mx3+w2*my1*w3*mx4-w2*my1*w4*mx3;
  return CGAL_NTS sign(coeff2);
}

template <typename NT>
inline
int coeff1_dor(const NT& mx1, const NT& nx1,
               const NT& my1, const NT& ny1, const NT& w1,
               const NT& mx2, const NT& nx2,
               const NT& my2, const NT& ny2, const NT& w2,
               const NT& mx3, const NT& nx3,
               const NT& my3, const NT& ny3, const NT& w3,
               const NT& mx4, const NT& nx4,
               const NT& my4, const NT& ny4, const NT& w4)
{
  NT coeff1 = -w1*my2*w3*nx4+w1*mx2*w3*ny4+w1*my2*w4*nx3-w1*mx2*w4*ny3+
               w1*nx2*w3*my4-w1*nx2*w4*my3+w2*my1*w3*nx4-w2*mx1*w3*ny4-
               w2*my1*w4*nx3+w2*mx1*w4*ny3-w2*nx1*w3*my4+w2*nx1*w4*my3-
               w1*ny2*w3*mx4+w1*ny2*w4*mx3+w2*ny1*w3*mx4-w2*ny1*w4*mx3;
  return CGAL_NTS sign(coeff1);
}

template <typename NT>
inline
int coeff0_dor(const NT& /*mx1*/, const NT& nx1,
               const NT& /*my1*/, const NT& ny1, const NT& w1,
               const NT& /*mx2*/, const NT& nx2,
               const NT& /*my2*/, const NT& ny2, const NT& w2,
               const NT& /*mx3*/, const NT& nx3,
               const NT& /*my3*/, const NT& ny3, const NT& w3,
               const NT& /*mx4*/, const NT& nx4,
               const NT& /*my4*/, const NT& ny4, const NT& w4)
{
  NT coeff0 = w1*nx2*w3*ny4-w1*nx2*w4*ny3-w2*nx1*w3*ny4+w2*nx1*w4*ny3-
              w1*ny2*w3*nx4+w1*ny2*w4*nx3+w2*ny1*w3*nx4-w2*ny1*w4*nx3;
  return CGAL_NTS sign(coeff0);
}

DEFCOUNTER(ord2)
DEFCOUNTER(ord1)
DEFCOUNTER(ord0)


template <typename RT>
inline
int orientation(const Extended_direction<RT>& d1,
                const Extended_direction<RT>& d2)
{
  Extended_point<RT> p1(d1.p1()), p2(d1.p2()),
                     p3(d2.p1()), p4(d2.p2());
  int res;
  try { INCTOTAL(ord2); Protect_FPU_rounding<true> Protection;
    res = coeff2_dor(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                     p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                     p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                     p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  } catch (Uncertain_conversion_exception&) { INCEXCEPTION(ord2);
    res = coeff2_dor(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                     p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                     p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                     p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  if ( res != 0 ) return res;

  try { INCTOTAL(ord1); Protect_FPU_rounding<true> Protection;
    res = coeff1_dor(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                     p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                     p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                     p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  } catch (Uncertain_conversion_exception&) { INCEXCEPTION(ord1);
    res = coeff1_dor(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                     p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                     p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                     p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  if ( res != 0 ) return res;
  try { INCTOTAL(ord0); Protect_FPU_rounding<true> Protection;
    res = coeff0_dor(p1.mxD(),p1.nxD(),p1.myD(),p1.nyD(),p1.hwD(),
                     p2.mxD(),p2.nxD(),p2.myD(),p2.nyD(),p2.hwD(),
                     p3.mxD(),p3.nxD(),p3.myD(),p3.nyD(),p3.hwD(),
                     p4.mxD(),p4.nxD(),p4.myD(),p4.nyD(),p4.hwD());
  } catch (Uncertain_conversion_exception&) { INCEXCEPTION(ord0);
    res = coeff0_dor(p1.mx(),p1.nx(),p1.my(),p1.ny(),p1.hw(),
                     p2.mx(),p2.nx(),p2.my(),p2.ny(),p2.hw(),
                     p3.mx(),p3.nx(),p3.my(),p3.ny(),p3.hw(),
                     p4.mx(),p4.nx(),p4.my(),p4.ny(),p4.hw());
  }
  return res;
}

template <typename RT>
inline
bool operator==(const Extended_direction<RT>& d1,
                const Extended_direction<RT>& d2)
{
  return orientation(d1,d2) == 0 &&
         d1.dx_sign() == d2.dx_sign() &&
         d1.dy_sign() == d2.dy_sign();
}

template <typename RT>
inline
bool operator!=(const Extended_direction<RT>& d1,
                const Extended_direction<RT>& d2)
{ return !(d1==d2); }


template <typename RT>
bool strictly_ordered_ccw(const Extended_direction<RT>& d1,
                          const Extended_direction<RT>& d2,
                          const Extended_direction<RT>& d3)
{
  if (d1 == d3) return (d1 != d2);
  int or12 = orientation(d1,d2);
  int or13 = orientation(d1,d3);
  int or32 = orientation(d3,d2);
  if ( or13 >= 0 ) // not right_turn
    return ( or12 > 0 && or32 < 0 );
  else // ( or13 < 0 ) right_turn
    return ( or12 > 0 || or32 < 0 );
}

template <typename RT>
inline
bool operator<(const Extended_direction<RT>& d1,
               const Extended_direction<RT>& d2)
{ Extended_direction<RT> d0(1,0);
  bool d0d1eq = (d1 == d0);
  bool d0d2eq = (d2 == d0);
  return ( (d0d1eq && !d0d2eq) ||
           ( strictly_ordered_ccw(d0,d1,d2) && (! d0d2eq) ) );
}



template<class Kernel>
struct Is_extended_kernel;
template <typename RT_>
class Filtered_extended_homogeneous;

template<class T>
struct Is_extended_kernel<Filtered_extended_homogeneous<T> > {
       typedef Tag_true value_type;
};


template <typename RT_>
class Filtered_extended_homogeneous {
typedef Filtered_extended_homogeneous<RT_> Self;

public:
typedef CGAL::Homogeneous<RT_> Standard_kernel;
typedef typename Standard_kernel::RT           Standard_RT;
typedef typename Standard_kernel::FT           Standard_FT;
typedef typename Standard_kernel::Point_2      Standard_point_2;
typedef typename Standard_kernel::Segment_2    Standard_segment_2;
typedef typename Standard_kernel::Line_2       Standard_line_2;
typedef typename Standard_kernel::Direction_2  Standard_direction_2;
typedef typename Standard_kernel::Ray_2        Standard_ray_2;
typedef typename Standard_kernel::Aff_transformation_2
  Standard_aff_transformation_2;

typedef SPolynomial<RT_>        RT;
typedef SQuotient<RT_>          FT;
typedef Extended_point<RT_>     Point_2;
typedef Extended_segment<RT_>   Segment_2;
typedef Extended_direction<RT_> Direction_2;
#ifdef KERNEL_CHECK
typedef Extended_homogeneous<RT_>         CheckKernel;
typedef typename CheckKernel::Point_2     CheckPoint;
typedef typename CheckKernel::Direction_2 CheckDirection;
typedef typename CheckKernel::Segment_2   CheckSegment;
CheckKernel K;

CheckSegment convert(const Segment_2& s) const
{ return CheckSegment(s.source().checkrep(),
                      s.target().checkrep()); }
CheckDirection convert(const Direction_2& d) const
{ return K.construct_direction(d.p2().checkrep(),d.p1().checkrep()); }

#endif // KERNEL_CHECK




enum Point_type { SWCORNER=1, LEFTFRAME, NWCORNER,
                  BOTTOMFRAME, STANDARD, TOPFRAME,
                  SECORNER, RIGHTFRAME, NECORNER };

Standard_RT dx(const Standard_line_2& l) const { return l.b(); }
Standard_RT dy(const Standard_line_2& l) const { return -l.a(); }
Standard_FT abscissa_distance(const Standard_line_2& l) const {
  typename CGAL::Rational_traits<typename Standard_kernel::FT> rat_traits;
  return rat_traits.make_rational(-l.c(), l.b());
  //  return -l.c() / l.b();
}

Point_type determine_type(const Standard_line_2& l) const
{
  // CGAL_NEF_TRACEN("determine_type "<<l);
  Standard_RT adx = CGAL_NTS abs(dx(l)), ady = CGAL_NTS abs(dy(l));
  int sdx = CGAL_NTS sign(dx(l)), sdy = CGAL_NTS sign(dy(l));
  int cmp_dx_dy = CGAL_NTS compare(adx,ady), s(1);
  // CGAL_NEF_TRACEN("   "<<cmp_dx_dy<<" "<<sdx<<" "<<sdy);
  if (sdx < 0 && ( cmp_dx_dy > 0 || ( cmp_dx_dy == 0 &&
                                      sdy != (s=CGAL_NTS sign(abscissa_distance(l)))))) {
    if (0 == s) return ( sdy < 0 ? SWCORNER : NWCORNER );
    else        return LEFTFRAME;
  } else if (sdx > 0 && ( cmp_dx_dy > 0 || ( cmp_dx_dy == 0 &&
                                             sdy != (s=CGAL_NTS sign(abscissa_distance(l)))))) {
    if (0 == s) return ( sdy < 0 ? SECORNER : NECORNER );
    else        return RIGHTFRAME;
  } else if (sdy < 0 && ( cmp_dx_dy < 0 || ( cmp_dx_dy == 0 &&
                                             abscissa_distance(l) < Standard_FT(0)))) {
    return BOTTOMFRAME;
  } else if (sdy > 0 && ( cmp_dx_dy < 0 || ( cmp_dx_dy == 0 &&
                                             abscissa_distance(l) > Standard_FT(0)))) {
    return TOPFRAME;
  }
  CGAL_error_msg(" determine_type: degenerate line.");
  return (Point_type)-1; // never come here
}

Point_2 epoint(const Standard_RT& m1, const Standard_RT& n1,
                  const Standard_RT& m2, const Standard_RT& n2,
                                 const Standard_RT& n3) const
{ return Point_2(m1,n1,m2,n2,n3); }

public:

Point_2 construct_point(const Standard_point_2& p) const
{ return Point_2(p.hx(), p.hy(), p.hw()); }

Point_2 construct_point(const Standard_line_2& l, Point_type& t) const
{
  t = determine_type(l);
  // CGAL_NEF_TRACEN("construct_point(line)"<<l<<" "<<t);
  Point_2 res;
  switch (t) {
    case SWCORNER:   res = epoint(-1, 0, -1, 0, 1); break;
    case NWCORNER:   res = epoint(-1, 0,  1, 0, 1); break;
    case SECORNER:   res = epoint( 1, 0, -1, 0, 1); break;
    case NECORNER:   res = epoint( 1, 0,  1, 0, 1); break;
    case LEFTFRAME:  res = epoint(-l.b(), 0,  l.a(), -l.c(), l.b());
                     break;
    case RIGHTFRAME: res = epoint( l.b(), 0, -l.a(), -l.c(), l.b());
                     break;
    case BOTTOMFRAME: res = epoint( l.b(), -l.c(), -l.a(), 0, l.a());
                     break;
    case TOPFRAME: res = epoint(-l.b(), -l.c(),  l.a(), 0, l.a());
                     break;
    default: CGAL_error_msg("EPoint type not correct!");
  }
  return res;
}

Point_2 construct_point(const Standard_point_2& p1,
                        const Standard_point_2& p2,
                        Point_type& t) const
{ return construct_point(Standard_line_2(p1,p2),t); }
Point_2 construct_point(const Standard_line_2& l) const
{ Point_type dummy; return construct_point(l,dummy); }
Point_2 construct_point(const Standard_point_2& p1,
                        const Standard_point_2& p2) const
{ return construct_point(Standard_line_2(p1,p2)); }
Point_2 construct_point(const Standard_point_2& p,
                        const Standard_direction_2& d) const
{ return construct_point(Standard_line_2(p,d)); }
Point_2 construct_opposite_point(const Standard_line_2& l) const
{ Point_type dummy; return construct_point(l.opposite(),dummy); }

Point_type type(const Point_2& p) const
{
  if (p.is_standard()) return STANDARD;
  // now we are on the square frame
  RT rx = p.hx();
  RT ry = p.hy();
  int sx = CGAL_NTS sign(rx);
  int sy = CGAL_NTS sign(ry);
  if (sx < 0) rx = -rx;
  if (sy < 0) ry = -ry;
  if (rx>ry) {
    if (sx > 0) return RIGHTFRAME;
    else        return LEFTFRAME;
  }
  if (rx<ry) {
    if (sy > 0) return TOPFRAME;
    else        return BOTTOMFRAME;
  }
  // now (rx == ry)
  if (sx==sy) {
    if (sx < 0) return SWCORNER;
    else        return NECORNER;
  } else { CGAL_assertion(sx==-sy);
    if (sx < 0) return NWCORNER;
    else        return SECORNER;
  }
}


bool is_standard(const Point_2& p) const
{ return p.is_standard();  }

Standard_point_2 standard_point(const Point_2& p) const
{ CGAL_assertion(is_standard(p));
  return Standard_point_2(p.nx(),p.ny(),p.hw());
}

Standard_line_2 standard_line(const Point_2& p) const
{ CGAL_assertion(!p.is_standard());
  Standard_point_2 p0(p.nx(),p.ny(),p.hw());
  Standard_point_2 p1(p.mx()+p.nx(),p.my()+p.ny(),p.hw());
  return Standard_line_2(p0,p1);
}

Standard_ray_2 standard_ray(const Point_2& p) const
{ CGAL_assertion(!p.is_standard());
  Standard_line_2 l = standard_line(p);
  Standard_direction_2 d = l.direction();
  Standard_point_2 q = l.point(0);
  return Standard_ray_2(q,d);
}

Point_2 NE() const { return construct_point(Standard_line_2(-1, 1,0)); }
Point_2 SE() const { return construct_point(Standard_line_2( 1, 1,0)); }
Point_2 NW() const { return construct_point(Standard_line_2(-1,-1,0)); }
Point_2 SW() const { return construct_point(Standard_line_2( 1,-1,0)); }

int orientation(const Point_2& p1, const Point_2& p2, const Point_2& p3)
const
{ CHECK(K.orientation(p1.checkrep(),p2.checkrep(),p3.checkrep()),
        CGAL::orientation(p1,p2,p3))
  return CGAL::orientation(p1,p2,p3); }

bool left_turn(const Point_2& p1, const Point_2& p2, const Point_2& p3)
const
{ return orientation(p1,p2,p3) > 0; }

bool first_pair_closer_than_second(
  const Point_2& p1, const Point_2& p2,
  const Point_2& p3, const Point_2& p4) const
{ CHECK(K.first_pair_closer_than_second(p1.checkrep(),p2.checkrep(),
                                        p3.checkrep(),p4.checkrep()),
        CGAL::compare_pair_dist(p1,p2,p3,p4)<0)
  return CGAL::compare_pair_dist(p1,p2,p3,p4)<0; }

int compare_xy(const Point_2& p1, const Point_2& p2) const
{ CHECK(K.compare_xy(p1.checkrep(),p2.checkrep()),
        CGAL::compare_xy(p1,p2))
  return CGAL::compare_xy(p1,p2); }

int compare_x(const Point_2& p1, const Point_2& p2) const
{ CHECK(K.compare_x(p1.checkrep(),p2.checkrep()),
        CGAL::compare_x(p1,p2))
  return CGAL::compare_x(p1,p2); }

int compare_y(const Point_2& p1, const Point_2& p2) const
{ CHECK(K.compare_y(p1.checkrep(),p2.checkrep()),
        CGAL::compare_y(p1,p2))
  return CGAL::compare_y(p1,p2); }

bool strictly_ordered_along_line(
  const Point_2& p1, const Point_2& p2, const Point_2& p3) const
{ CHECK(K.strictly_ordered_along_line(
          p1.checkrep(),p2.checkrep(),p3.checkrep()),
        CGAL::strictly_ordered_along_line(p1,p2,p3))
  return CGAL::strictly_ordered_along_line(p1,p2,p3); }


Segment_2 construct_segment(const Point_2& p, const Point_2& q) const
{ return Segment_2(p,q); }

Point_2 source(const Segment_2& s) const
{ return s.source(); }

Point_2 target(const Segment_2& s) const
{ return s.target(); }

bool is_degenerate(const Segment_2& s) const
{ return s.source()==s.target(); }

int orientation(const Segment_2& s, const Point_2& p) const
{ return orientation(s.source(),s.target(),p); }

Point_2 intersection(const Segment_2& s1, const Segment_2& s2) const
{ CHECK(CGAL::intersection(s1,s2).checkrep(),
        K.intersection(convert(s1),convert(s2)))
  return CGAL::intersection(s1,s2); }

bool contains(const Segment_2& s, const Point_2& p) const
/*{\Mop returns true iff |s| contains |p|.}*/
{ return CGAL::contains(s,p); }

Direction_2 construct_direction(
  const Point_2& p1, const Point_2& p2) const
{ return Direction_2(p1,p2); }

bool strictly_ordered_ccw(const Direction_2& d1,
  const Direction_2& d2, const Direction_2& d3) const
{ CHECK(K.strictly_ordered_ccw(convert(d1),convert(d2),convert(d3)),
        CGAL::strictly_ordered_ccw(d1,d2,d3));
  return CGAL::strictly_ordered_ccw(d1,d2,d3); }

void print_statistics() const
{
  std::cout << "Statistics of filtered kernel:\n";
  std::cout << "total failed double filter stages = (now needs CGAL_PROFILE)\n";
  PRINT_CHECK_ENABLED;
  PRINT_STATISTICS(or2);
  PRINT_STATISTICS(or1);
  PRINT_STATISTICS(or0);
  PRINT_STATISTICS(cmpx1);
  PRINT_STATISTICS(cmpx0);
  PRINT_STATISTICS(cmpy1);
  PRINT_STATISTICS(cmpy0);
  PRINT_STATISTICS(cmppd2);
  PRINT_STATISTICS(cmppd1);
  PRINT_STATISTICS(cmppd0);
  PRINT_STATISTICS(ord2);
  PRINT_STATISTICS(ord1);
  PRINT_STATISTICS(ord0);
}

template <class Forward_iterator>
void determine_frame_radius(Forward_iterator start, Forward_iterator end,
                            Standard_RT& R0) const
{ Standard_RT R;
  while ( start != end ) {
    Point_2 p = *start;
    if ( is_standard(p) ) {
      R = (CGAL::max)(CGAL_NTS abs(p.mx())/p.hw(),
                      CGAL_NTS abs(p.my())/p.hw());
    } else {
      RT rx = CGAL_NTS abs(p.hx()), ry = CGAL_NTS abs(p.hy());
      if ( rx[1] > ry[1] )      R = CGAL_NTS abs(ry[0]-rx[0])/(rx[1]-ry[1]);
      else if ( rx[1] < ry[1] ) R = CGAL_NTS abs(rx[0]-ry[0])/(ry[1]-rx[1]);
      else /* rx[1] == ry[1] */ R = CGAL_NTS abs(rx[0]-ry[0])/(2*p.hw());
    }
    R0 = (CGAL::max)(R+1,R0); ++start;
  }
}

const char* output_identifier() const
{ return "Filtered_extended_homogeneous"; }



};



} //namespace CGAL

#undef CHECK
#undef KERNEL_CHECK
#undef REDUCE_INTERSECTION_POINTS
#undef KERNEL_ANALYSIS
#undef COUNTER
#undef INCTOTAL
#undef INCEXCEPTION
#undef PRINT_STATISTICS
#undef PRINT_CHECK_ENABLED

#include <CGAL/enable_warnings.h>

#endif // CGAL_FILTERED_EXTENDED_HOMOGENEOUS_H
