#ifndef INPUT_TRAITS_H
#define INPUT_TRAITS_H

template <class My_NT>
struct Input_traits {};

/*! int */
template <>
struct Input_traits<int> {
  typedef int                                   Input_int_type;
  typedef float                                 Input_float_type;
  typedef int                                   Input_rat_type;
};

/*! float */
template <>
struct Input_traits<float> {
  typedef int                                   Input_int_type;
  typedef float                                 Input_float_type;
  typedef float                                 Input_rat_type;
};

/*! double */
template <>
struct Input_traits<double> {
  typedef int                                   Input_int_type;
  typedef double                                Input_float_type;
  typedef double                                Input_rat_type;
};

#ifdef CGAL_MP_FLOAT_H
/*! CGAL::MP_Float */
template <>
struct Input_traits<CGAL::MP_Float> {
  typedef CGAL::MP_Float                        Input_int_type;
  typedef CGAL::MP_Float                        Input_float_type;
  typedef CGAL::MP_Float                        Input_rat_type;
};
#endif

#ifdef __GMP_PLUSPLUS__
/*! ::mpz_class */
template <>
struct Input_traits<::mpz_class> {
  typedef ::mpz_class                           Input_int_type;
  typedef ::mpz_class                           Input_float_type;
  typedef ::mpz_class                           Input_rat_type;
};
#endif

/*! CGAL::Gmpz */
template <>
struct Input_traits<CGAL::Gmpz> {
  typedef CGAL::Gmpz                            Input_int_type;
  typedef CGAL::Gmpz                            Input_float_type;
  typedef CGAL::Gmpz                            Input_rat_type;
};

/*! leda_integer */
template <>
struct Input_traits<leda_integer> {
  typedef leda_integer                          Input_int_type;
  typedef leda_integer                          Input_float_type;
  typedef leda_integer                          Input_rat_type;
};

#ifdef _CORE_BIGINT_H_
/*! CORE::BigInt */
template <>
struct Input_traits<CORE::BigInt> {
  typedef CORE::BigInt                          Input_int_type;
  typedef CORE::BigInt                          Input_float_type;
  typedef CORE::BigInt                          Input_rat_type;
};
#endif

#ifdef CGAL_MP_FLOAT_H
/*! CGAL::Quotient<CGAL::MP_Float> */
template <>
struct Input_traits<CGAL::Quotient<CGAL::MP_Float> > {
  typedef CGAL::MP_Float                        Input_int_type;
  typedef CGAL::MP_Float                        Input_float_type;
  typedef CGAL::Quotient<CGAL::MP_Float>        Input_rat_type;
};
#endif

#ifdef __GMP_PLUSPLUS__
/*! CGAL::Quotient<::mpz_struct> */
template <>
struct Input_traits<CGAL::Quotient<::mpz_struct> > {
  typedef ::mpz_struct                          Input_int_type;
  typedef ::mpz_struct                          Input_float_type;
  typedef ::mpz_struct                          Input_rat_type;
};
#endif

/*! CGAL::Quotient<CGAL::Gmpz> */
template <>
struct Input_traits<CGAL::Quotient<CGAL::Gmpz> > {
  typedef CGAL::Gmpz                            Input_int_type;
  typedef CGAL::Gmpz                            Input_float_type;
  typedef CGAL::Gmpz                            Input_rat_type;
};

#ifdef __GMP_PLUSPLUS__
/*! ::mpq_struct */
template <>
struct Input_traits<::mpq_struct> {
  typedef ::mpz_struct                          Input_int_type;
  typedef ::mpz_struct                          Input_float_type;
  typedef ::mpz_struct                          Input_rat_type;
};
#endif

/*! CGAL::Gmpq */
template <>
struct Input_traits<CGAL::Gmpq> {
  typedef CGAL::Gmpz                            Input_int_type;
  typedef CGAL::Gmpz                            Input_float_type;
  typedef CGAL::Gmpz                            Input_rat_type;
};

/*! leda_rational */
template <>
struct Input_traits<leda_rational> {
  typedef leda_integer                          Input_int_type;
  typedef leda_integer                          Input_float_type;
  typedef leda_integer                          Input_rat_type;
};

#ifdef CGAL_LAZY_EXACT_NT_H
/*! CGAL::Lazy_exact_nt<leda_rational> */
template <>
struct Input_traits<CGAL::Lazy_exact_nt<leda_rational> > {
  typedef leda_integer                          Input_int_type;
  typedef leda_integer                          Input_float_type;
  typedef leda_integer                          Input_rat_type;
};

/*! CGAL::Lazy_exact_nt<leda_rational> */
template <>
struct Input_traits<CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > > {
  typedef CGAL::MP_Float                        Input_int_type;
  typedef CGAL::MP_Float                        Input_float_type;
  typedef CGAL::Quotient<CGAL::MP_Float>        Input_rat_type;

  typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  Output_type;
};
#endif

#endif
