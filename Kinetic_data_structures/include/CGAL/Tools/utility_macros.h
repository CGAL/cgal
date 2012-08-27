// Copyright (c) 2007 (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>


#ifndef CGAL_UTILITY_MACROS_H
#define CGAL_UTILITY_MACROS_H

#define CGAL_SUBSCRIPT(type, expr) type& operator[](unsigned int i){ expr;} \
  const type& operator[](unsigned int i) const { expr;}

#define CGAL_COPY_CONSTRUCTOR(TC) TC(const TC &o){copy_from(o);}\
  TC& operator=(const TC &o) {copy_from(o); return *this;}

// 
#define CGAL_GET(type, name, expr) const type & name() const{expr;}
#define CGAL_GETOBJECT(UC, lc, expr) const UC& lc##_object() const {expr;}
#define CGAL_GETNR(type, name, expr) type name() const{expr;}
//#define CGAL_GET(type, name, expr) typename boost::call_traits<type>::param_type name() const {expr;}

#define CGAL_IS(name, expr) bool is_##name() const {expr;}
#define CGAL_SET_IS(name, expr) void set_is_##name(bool tf) {expr;}

#define CGAL_SET(type, name, expr) void set_##name(const type &k) {expr;}

#define CGAL_GETSET(type, name, var)		\
  CGAL_GET(type, name, return var)              \
  CGAL_SET(type, name, var = k)

#define CGAL_OUTPUT(type)						\
  inline std::ostream& operator<<(std::ostream&out, const type &t){	\
    return t.write(out);						\
  }                                                                     \
  CGAL_REQUIRE_SEMICOLON_NAMESPACE

#define CGAL_OUTPUT1(type)						\
  template <class A>							\
  inline std::ostream& operator<<(std::ostream&out, const type<A> &t){	\
    return t.write(out);						\
  }

#define CGAL_OUTPUT2(T)							\
  template <class A, class B>						\
  inline std::ostream& operator<<(std::ostream&out, const T<A,B> &t){	\
    return t.write(out);						\
  }

#define CGAL_ITERATOR(uc_name, lc_name, const_it_type,it_type, bexpr, eexpr) \
  typedef boost::iterator_range< it_type > uc_name##s;                  \
  uc_name##s lc_name##s() {                                             \
    return uc_name##s(it_type(bexpr), it_type(eexpr));                  \
  }                                                                     \
  typedef boost::iterator_range< const_it_type > uc_name##_consts;      \
  uc_name##_consts lc_name##s() const {                                 \
    return uc_name##_consts(const_it_type(bexpr), const_it_type(eexpr)); \
  }

#define CGAL_CONST_ITERATOR(uc_name, lc_name, it_type, bexpr, eexpr)	\
  typedef boost::iterator_range< it_type > uc_name##s;                  \
  uc_name##s lc_name##s() const {                                       \
    return uc_name##s(bexpr, eexpr);                                    \
  }



#define CGAL_CONST_FIND(ucname, fexpr, expr)                       \
  bool contains(ucname##_key k) const {                            \
    return fexpr != eexpr;                                         \
  }                                                                \
  std::iterator_traits<ucname##s::iterator>::value_type get(ucname##_key k) const {      \
    CGAL_assertion(contains(k));                                   \
    return *fexpr;                                                 \
  }                                                                \
  ucname##s::const_iterator find(ucname##_key k) {                 \
    return fexpr;                                                  \
  }

#define CGAL_FIND(ucname, fexpr, eexpr)                         \
  bool contains(ucname##_key k) const {                         \
    return fexpr != eexpr;                                      \
  }                                                             \
  std::iterator_traits< ucname##s::iterator >::reference get(ucname##_key k) { \
    CGAL_assertion(contains(k));                                \
    return *fexpr;                                              \
  }                                                             \
  ucname##s::iterator find(ucname##_key k) {                    \
    return fexpr;                                               \
  }                                                             \
  std::iterator_traits< ucname##s::const_iterator >::value_type get(ucname##_key k) const {      \
    CGAL_assertion(contains(k));                                   \
    return *fexpr;                                                 \
  }                                                                \
  ucname##_consts::const_iterator find(ucname##_key k) const {           \
    return fexpr;                                                  \
  }



#define CGAL_INSERT(ucname, expr)					\
  void insert(ucname##_key k, const ucname &m) {expr;}

#define CGAL_INSERTNK(ucname, expr)			\
  void insert(const ucname &m) {expr;}

#define CGAL_SWAP(type)				\
  inline void swap(type &a, type &b) {		\
    a.swap_with(b);				\
  }

#define CGAL_SWAP1(type)			\
  template <class A>				\
 inline void swap(type<A> &a, type<A> &b) {	\
    a.swap_with(b);				\
  }

#define CGAL_SWAP2(type)			\
template <class A, class B>			\
 inline void swap(type<A,B> &a, type<A,B> &b) {	\
    a.swap_with(b);				\
}

#define CGAL_ISWAP(name)			\
  std::swap(name, o.name)


#define CGAL_IFNONEQUAL(a,b,cmp) if (a cmp b) return true;	\
  else if (b cmp a) return false

#define CGAL_COMPARISONS bool operator==(const This &o) const {		\
    return compare(o) == CGAL::EQUAL;					\
  }									\
  bool operator!=(const This &o) const {				\
    return compare(o) != CGAL::EQUAL;					\
  }									\
  bool operator<(const This &o) const {					\
    return compare(o) == CGAL::SMALLER;					\
  }									\
  bool operator>(const This &o) const {					\
    return compare(o) == CGAL::LARGER;					\
  }									\
  bool operator>=(const This &o) const {				\
    return compare(o) != CGAL::SMALLER;					\
  }									\
  bool operator<=(const This &o) const {				\
    return compare(o) != CGAL::LARGER;					\
  }

#define CGAL_COMPARISONS_COMPARE bool operator==(const This &o) const { \
    return compare(o)==0;						\
  }									\
  bool operator!=(const This &o) const {				\
    return compare(o)!=0;						\
  }									\
  bool operator<(const This &o) const {					\
    return compare(o) < 0;						\
  }									\
  bool operator>(const This &o) const {					\
    return compare(o) > 0;						\
  }									\
  bool operator>=(const This &o) const {				\
    return compare(o) >=0;						\
  }									\
  bool operator<=(const This &o) const {				\
    return compare(o) <= 0;                                             \
  }



#define CGAL_COMPARISONS1(field) bool operator==(const This &o) const { \
    return (field== o.field);						\
  }									\
  bool operator!=(const This &o) const {				\
    return (field!= o.field);						\
  }									\
  bool operator<(const This &o) const {					\
    return (field< o.field);						\
  }									\
  bool operator>(const This &o) const {					\
    return (field> o.field);						\
  }									\
  bool operator>=(const This &o) const {				\
    return (field>= o.field);						\
  }									\
  bool operator<=(const This &o) const {				\
    return (field<= o.field);						\
  }




#define CGAL_COMPARISONS2(a, b) bool operator==(const This &o) const {	\
    return (a== o.a && b== o.b);					\
  }									\
  bool operator!=(const This &o) const {				\
    return (a!= o.a || b != o.b);					\
  }									\
  bool operator<(const This &o) const {					\
    if (a< o.a ) return true;						\
    else if (a > o.a) return false;					\
    else return b < o.b;						\
  }									\
  bool operator>(const This &o) const {					\
    if (a> o.a ) return true;						\
    else if (a < o.a) return false;					\
    else return b > o.b;						\
  }									\
  bool operator>=(const This &o) const {				\
    return !operator<(o);						\
  }									\
  bool operator<=(const This &o) const {				\
    return !operator>(o);						\
  }


#define CGAL_COMPARISONS3(a, b, c)					\
  bool operator==(const This &o) const {				\
    return (a== o.a && b== o.b && c == o.c);				\
  }									\
  bool operator!=(const This &o) const {				\
    return (a!= o.a || b != o.b || c != o.c);				\
  }									\
  bool operator<(const This &o) const {					\
    if (a < o.a ) return true;						\
    else if (a > o.a) return false;					\
    else if (b < o.b) return true;					\
    else if (b > o.b) return false;					\
    else return c < o.c;						\
  }									\
  bool operator>(const This &o) const {					\
    if (a > o.a ) return true;						\
    else if (a < o.a) return false;					\
    else if (b > o.b) return true;					\
    else if (b < o.b) return false;					\
    else return c > o.c;						\
  }									\
  bool operator>=(const This &o) const {				\
    return !operator<(o);						\
  }									\
  bool operator<=(const This &o) const {				\
    return !operator>(o);						\
  }


#define CGAL_REAL_EMBEDDABLE_BODY					\
  class Abs								\
    : public std::unary_function< Type, Type > {                        \
  public:								\
    Type operator()( const Type& x ) const {				\
      if (x < Type(0)) return -x;					\
      else return x;							\
    }									\
  };									\
  									\
  class Sgn								\
    : public std::unary_function< Type, ::CGAL::Sign > {                \
  public:								\
    ::CGAL::Sign operator()( const Type& x ) const {			\
      return static_cast<CGAL::Sign>(x.compare(0));			\
    }									\
  };									\
  									\
  class Compare								\
    : public std::binary_function< Type, Type,				\
			      Comparison_result > {			\
  public:								\
      Comparison_result operator()( const Type& x,			\
				    const Type& y ) const {		\
	return x.compare(y);						\
      }									\
    									\
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,		\
							 Comparison_result ); \
									\
  };									\
									\
  class To_double							\
    : public std::unary_function< Type, double > {                      \
  public:								\
    double operator()( const Type& x ) const {				\
      return x.approximation(.00000001);				\
    }									\
  };									\
									\
  class To_interval							\
    : public std::unary_function< Type, std::pair< double, double > > {	\
  public:								\
    std::pair<double, double> operator()( const Type& x ) const {	\
									\
      return x.approximating_interval(.00001);				\
    }									\
  }			

#define CGAL_HAS_INFINITY_BODY						\
  static const bool is_specialized = true;				\
  static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw ()             \
  {return -T::infinity();}                                              \
  static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw ()             \
  {return T::infinity();}                                               \
  static const int digits =0;						\
  static const int digits10 =0;						\
  static const bool is_signed = true;					\
  static const bool is_integer = false;					\
  static const bool is_exact = true;					\
  static const int radix =0;						\
  static T epsilon() throw(){return T(0);}				\
  static T round_error() throw(){return T(0);}				\
  static const int min_exponent=0;					\
  static const int min_exponent10=0;					\
  static const int max_exponent=0;					\
  static const int max_exponent10=0;					\
  static const bool has_infinity=true;					\
  static const bool has_quiet_NaN = true;				\
  static const bool has_signaling_NaN= false;				\
  static const float_denorm_style has_denorm= denorm_absent;		\
  static const bool has_denorm_loss = false;				\
  static T infinity() throw() {return T::infinity();}			\
  static T quiet_NaN() throw(){return T();}				\
  static T denorm_min() throw() {return T(0);}				\
  static const bool is_iec559=false;					\
  static const bool is_bounded =false;					\
  static const bool is_modulo= false;					\
  static const bool traps = false;					\
  static const bool tinyness_before =false;				\
  static const float_round_style round_stype = round_toward_zero;	\
    

#define CGAL_REAL_EMBEDDABLE1(name)					\
  namespace CGAL {							\
  template <class T>							\
  class Real_embeddable_traits< name<T> >				\
    : public INTERN_RET::Real_embeddable_traits_base<name<T>, Tag_true>{ \
  public:								\
    typedef name<T>  Type;						\
    CGAL_REAL_EMBEDDABLE_BODY;						\
  };									\
  } //namespace CGAL							


#define CGAL_REAL_EMBEDDABLE2(name)					\
  namespace CGAL {							\
  template <class T, class A>						\
  class Real_embeddable_traits< name<T, A> >				\
    : public INTERN_RET::Real_embeddable_traits_base< name<T, A>, Tag_true>{ \
  public:								\
    typedef name<T, A>  Type;						\
    CGAL_REAL_EMBEDDABLE_BODY;						\
  };									\
  } //namespace CGAL							


#define CGAL_HAS_INFINITY1(name)					\
  namespace std								\
  {									\
    template <class Tr>							\
      class numeric_limits<name<Tr> >					\
    {									\
    public:								\
      typedef name<Tr> T;						\
      CGAL_HAS_INFINITY_BODY;						\
    };									\
  }




#define CGAL_REAL_EMBEDABLE2(name)					\
  namespace CGAL {							\
  template <class T,class U>						\
  class Real_embeddable_traits< name<T, U> >				\
    : public INTERN_RET::Real_embeddable_traits_base< name<T, U> , Tag_true>{ \
  public:								\
    typedef name<T, U>  Type;						\
    CGAL_REAL_EMBEDDABLE_BODY						\
  };									\
  } //namespace CGAL							

#define CGAL_HAS_INFINITY2(name)					\
  namespace std								\
  {									\
    template <class Tr, class Ur>					\
      class numeric_limits<name<Tr, Ur> >				\
    {									\
    public:								\
      typedef name<Tr, Ur> T;						\
      CGAL_HAS_INFINITY_BODY;						\
    };									\
  }

#endif
