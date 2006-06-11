#ifndef NUMBER_TYPE_TRAITS_H
#define NUMBER_TYPE_TRAITS_H

#include "number_type.h"

template <class My_NT>
struct Number_type_traits {};

/*! int */
template <>
struct Number_type_traits<int> {
  typedef int                                   RT;
  typedef int                                   FT;
};

/*! float */
template <>
struct Number_type_traits<float> {
  typedef float                                 RT;
  typedef float                                 FT;
};

#ifdef DOUBLE_H
/*! double */
template <>
struct Number_type_traits<CGAL::Double> {
  typedef double                                RT;
  typedef double                                FT;
};
#endif

#ifdef CGAL_MP_FLOAT_H
/*! CGAL::MP_Float */
template <>
struct Number_type_traits<CGAL::MP_Float> {
  typedef CGAL::MP_Float                        RT;
  typedef CGAL::MP_Float                        FT;
};
#endif

#if defined(__GMP_PLUSPLUS__) && 0
/*! ::mpz_class */
template <>
struct Number_type_traits<::mpz_class> {
  typedef ::mpz_class                           RT;
  typedef ::mpz_class                           FT;
};
#endif

#ifdef CGAL_GMPZ_H
/*! CGAL::Gmpz */
template <>
struct Number_type_traits<CGAL::Gmpz> {
  typedef CGAL::Gmpz                            RT;
  typedef CGAL::Gmpz                            FT;
};
#endif

#ifdef CGAL_LEDA_INTEGER_H
/*! leda_integer */
template <>
struct Number_type_traits<leda_integer> {
  typedef leda_integer                          RT;
  typedef leda_integer                          FT;
};
#endif

#ifdef _CORE_BIGINT_H_
/*! CORE::BigInt */
template <>
struct Number_type_traits<CORE::BigInt> {
  typedef CORE::BigInt                          RT;
  typedef CORE::BigInt                          FT;
};
#endif

#ifdef CGAL_MP_FLOAT_H
/*! CGAL::Quotient<CGAL::MP_Float> */
template <>
struct Number_type_traits<CGAL::Quotient<CGAL::MP_Float> > {
  typedef CGAL::MP_Float                        RT;
  typedef CGAL::Quotient<CGAL::MP_Float>        FT;
};
#endif

#if defined(__GMP_PLUSPLUS__) && 0
/*! CGAL::Quotient<::mpz_struct> */
template <>
struct Number_type_traits<CGAL::Quotient<::mpz_struct> > {
  typedef ::mpz_struct                          RT;
  typedef CGAL::Quotient<::mpz_struct>          FT;
};
#endif

#ifdef CGAL_GMPZ_H
/*! CGAL::Quotient<CGAL::Gmpz> */
template <>
struct Number_type_traits<CGAL::Quotient<CGAL::Gmpz> > {
  typedef CGAL::Gmpz                            RT;
  typedef CGAL::Quotient<CGAL::Gmpz>            FT;
};
#endif

#if defined(__GMP_PLUSPLUS__) && 0
/*! ::mpq_struct */
template <>
struct Number_type_traits<::mpq_struct> {
  typedef ::mpz_struct                          RT;
  typedef ::mpq_struct                          FT;
};
#endif

#if defined(CGAL_GMPQ_H)
/*! CGAL::Gmpq */
template <>
struct Number_type_traits<CGAL::Gmpq> {
  typedef CGAL::Gmpz                            RT;
  typedef CGAL::Gmpq                            FT;
};
#endif

#ifdef CGAL_LEDA_RATIONAL_H
/*! leda_rational */
template <>
struct Number_type_traits<leda_rational> {
  typedef leda_integer                          RT;
  typedef leda_rational                         FT;
};
#endif

#ifdef CGAL_LAZY_EXACT_NT_H
#ifdef CGAL_LEDA_RATIONAL_H
/*! CGAL::Lazy_exact_nt<leda_rational> */
template <>
struct Number_type_traits<CGAL::Lazy_exact_nt<leda_rational> > {
  typedef leda_integer                          RT;
  typedef leda_rational                         FT;
};
#endif

/*! CGAL::Lazy_exact_nt<CGAL::Gmpq> */
template <>
struct Number_type_traits<CGAL::Lazy_exact_nt<CGAL::Gmpq> > {
  typedef CGAL::Gmpz                            RT;
  typedef CGAL::Gmpq                            FT;
};

#ifdef CGAL_MP_FLOAT_H
/*! CGAL::Lazy_exact_nt<Quotient<MP_float>> */
template <>
struct Number_type_traits<CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > > {
  typedef CGAL::MP_Float                        RT;
  typedef CGAL::Quotient<CGAL::MP_Float>        FT;
};
#endif
#endif

//#ifdef _CORE_EXPR_H_  
#ifdef CGAL_CORE_EXPR_H
/*! CORE::Expr */
template <>
struct Number_type_traits<CORE::Expr> {
  typedef CORE::BigInt                          RT;
  typedef CORE::Expr                            FT;
};
#endif

#endif
