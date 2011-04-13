// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Quotient.C
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     :
//
//  The template class Quotient<NumberType> is based on the LEDA class
//  leda_rational written by Stefan Naeher and Christian Uhrig.
//  It is basically a templated version with restricted functionality
//  of the version of rational in LEDA release 3.3.
//  The modification was done by Stefan.Schirra@mpi-sb.mpg.de
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_QUOTIENT_C
#define CGAL_QUOTIENT_C

#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

template <class NumberType>
CGAL_KERNEL_MEDIUM_INLINE
Quotient<NumberType>&
Quotient<NumberType>::normalize()
{
  typedef NumberType NT;
  if (num == den)
  {
      num = den = NT(1);
      return (*this);
  }
  if (-num == den)
  {
      num = NT(-1);
      den = NT( 1);
      return (*this);
  }
  NumberType ggt = gcd(num, den);
  if (ggt != NT(1) )
  {
      num /= ggt;
      den /= ggt;
  }
  return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator+= (const Quotient<NumberType>& r)
{
    num = num * r.den + r.num * den;
    den *= r.den;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator-= (const Quotient<NumberType>& r)
{
    num = num * r.den - r.num * den;
    den *= r.den;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator*= (const Quotient<NumberType>& r)
{
    num *= r.num;
    den *= r.den;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator/= (const Quotient<NumberType>& r)
{
    typedef NumberType NT;
    CGAL_kernel_precondition( r.num != NT(0) );
    num *= r.den;
    den *= r.num;
    return (*this);
}
template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator+= (const NumberType& r)
{
    num = num + r * den;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator-= (const NumberType& r)
{
    num = num - r * den;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator*= (const NumberType& r)
{
    num *= r ;
    return (*this);
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>&
Quotient<NumberType>::operator/= (const NumberType& r)
{
    CGAL_kernel_precondition_code( typedef NumberType NT; )
    CGAL_kernel_precondition( r != NT(0) );
    den *= r ;
    return (*this);
}
template <class NumberType>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
quotient_cmp(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
    // In contrast to LEDA class rational, no assumptions
    // on the sign of  den  are made

    // code assumes that SMALLER == - 1;
    CGAL_kernel_precondition( SMALLER == static_cast<Comparison_result>(-1) );

    int xsign = CGAL_NTS sign(x.num) * CGAL_NTS sign(x.den) ;
    int ysign = CGAL_NTS sign(y.num) * CGAL_NTS sign(y.den) ;
    if (xsign == 0) return static_cast<Comparison_result>(-ysign);
    if (ysign == 0) return static_cast<Comparison_result>(xsign);
    // now (x != 0) && (y != 0)
    int diff = xsign - ysign;
    if (diff == 0)
    {
        int msign = CGAL_NTS sign(x.den) * CGAL_NTS sign(y.den);
        NumberType leftop  = x.num * y.den * NumberType(msign);
        NumberType rightop = y.num * x.den * NumberType(msign);
        if (leftop < rightop)
        {
            return SMALLER;
        }
        else
        {
            return (rightop < leftop) ? LARGER : EQUAL;
        }
    }
    else
    {
        return (xsign < ysign) ? SMALLER : LARGER;
    }
}

template <class NumberType>
std::ostream&
operator<<(std::ostream& s, const Quotient<NumberType>& r)
{
   s << r.numerator() << "/" << r.denominator();
   return s;
}

template <class NumberType>
std::istream&
operator>>(std::istream& in, Quotient<NumberType>& r)
{
  /* format  num/den  or simply  num  */

  char c = 0;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if ( !in ) return in;
  in.putback(c);

  NumberType num;
  NumberType den(1);
  in >> num;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if (( in ) && ( c == '/'))
  {
#ifndef CGAL_CFG_NO_LOCALE
      while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
      while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
      CGAL_kernel_assertion( in );
      in.putback(c);
      in >> den;
  }
  else
  {
      in.putback(c);
      if ( in.eof() ) in.clear();
  }
  r = Quotient<NumberType>( num, den);
  return in;
}

template <class NumberType>
inline
io_Operator
io_tag(const Quotient<NumberType>&)
{ return io_Operator(); }

template <class NumberType>
inline
NumberType
Quotient<NumberType>::numerator() const
{ return num; }

template <class NumberType>
inline
NumberType
Quotient<NumberType>::denominator() const
{ return den; }

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator+(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z = x;
  return z += y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator-(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{ return (Quotient<NumberType>(x) -= y); }

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator*(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z = x;
  return z *= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator/(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z = x;
  return z /= y;
}

template <class NumberType>
inline
Quotient<NumberType>
operator-(const Quotient<NumberType>& x)
{ return Quotient<NumberType>(-x.num,x.den); }

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator+(const NumberType& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z(x);
  return z += y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator+(const Quotient<NumberType>& x, const NumberType& y)
{
  Quotient<NumberType> z = x;
  return z += y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator-(const NumberType& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z(x);
  return z -= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator-(const Quotient<NumberType>& x, const NumberType& y)
{
  Quotient<NumberType> z = x;
  return z -= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator*(const NumberType& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z(x);
  return z *= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator*(const Quotient<NumberType>& x, const NumberType& y)
{
  Quotient<NumberType> z = x;
  return z *= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator/(const NumberType& x, const Quotient<NumberType>& y)
{
  Quotient<NumberType> z(x) ;
  return z /= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
Quotient<NumberType>
operator/(const Quotient<NumberType>& x, const NumberType& y)
{
  Quotient<NumberType> z = x;
  return z /= y;
}

template <class NumberType>
CGAL_KERNEL_INLINE
NumberType
quotient_truncation(const Quotient<NumberType>& r)
{ return (r.num / r.den); }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator==(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{ return x.num * y.den == x.den * y.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(x,y) == SMALLER; // < 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator==(const Quotient<NumberType>& x, const NumberType& y)
{ return x.den * y == x.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator==(const NumberType& x, const Quotient<NumberType>& y)
{ return y.den * x == y.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<(const Quotient<NumberType>& x, const NumberType& y)
{
  return quotient_cmp(x,Quotient<NumberType>(y)) == SMALLER; // < 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<(const NumberType& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(Quotient<NumberType>(x),y) == SMALLER; // < 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator!=(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{ return x.num * y.den != x.den * y.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<=(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(x,y) != LARGER; // <= 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator!=(const Quotient<NumberType>& x, const NumberType& y)
{ return x.den * y != x.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator!=(const NumberType& x, const Quotient<NumberType>& y)
{ return y.den * x != y.num; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<=(const Quotient<NumberType>& x, const NumberType& y)
{
  return quotient_cmp(x,Quotient<NumberType>(y)) != LARGER; // <= 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator<=(const NumberType& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(Quotient<NumberType>(x),y) != LARGER; // <= 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>(const Quotient<NumberType>& x, const NumberType& y)
{
  return  quotient_cmp(x,Quotient<NumberType>(y)) == LARGER; // > 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>(const NumberType& x, const Quotient<NumberType>& y)
{
  return  quotient_cmp(Quotient<NumberType>(x),y) == LARGER; // > 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>=(const Quotient<NumberType>& x, const NumberType& y)
{
  return quotient_cmp(x,Quotient<NumberType>(y)) != SMALLER; // >= 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>=(const NumberType& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(Quotient<NumberType>(x),y) != SMALLER; // >= 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  return  quotient_cmp(x,y) == LARGER; // > 0
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
operator>=(const Quotient<NumberType>& x, const Quotient<NumberType>& y)
{
  return quotient_cmp(x,y) != SMALLER; // >= 0
}

template <class NumberType>
double
to_double(const Quotient<NumberType>& q)   /* TODO */
{
  if (q.num == NumberType(0) )
  { return 0; }

  double nd = CGAL::to_double( q.num );

  if (q.den == NumberType(1) )
  { return nd; }

  double dd = CGAL::to_double( q.den );

  if ( is_finite( q.den )&&(is_finite( q.num )) )
  { return nd/dd ; }
  if ( CGAL_NTS abs(q.num) > CGAL_NTS abs(q.den) )
  {
      NumberType  nt_div = q.num / q.den;
      double divd = CGAL::to_double(nt_div);
      if ( divd >= ldexp(1.0,53) )
      { return divd; }
  }
  if ( CGAL_NTS abs(q.num) < CGAL_NTS abs(q.den) )
  { return 1.0 / CGAL::to_double( NumberType(1) / q ); }

  return nd/dd ;
}

template <class NumberType>
CGAL_KERNEL_INLINE
bool
is_valid(const Quotient<NumberType>& q)
{ return is_valid(q.num) && is_valid(q.den) ; }

template <class NumberType>
CGAL_KERNEL_INLINE
bool
is_finite(const Quotient<NumberType>& q)
{ return (is_finite(q.num))&&(is_finite(q.den)) ; }

template <class NumberType>
CGAL_KERNEL_INLINE
NumberType
denominator(const Quotient<NumberType>& q)
{ return q.den ; }

template <class NumberType>
CGAL_KERNEL_INLINE
NumberType
numerator(const Quotient<NumberType>& q)
{ return q.num ; }

template <class NumberType>
inline
Quotient_tag
number_type_tag(const Quotient<NumberType>& q)
{ return Quotient_tag(); }

template <class NumberType>
NumberType
gcd(const NumberType&, const NumberType&)
{ return NumberType(1); }

CGAL_END_NAMESPACE

#endif  // CGAL_QUOTIENT_C
