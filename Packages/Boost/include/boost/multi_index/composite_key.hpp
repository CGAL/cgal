/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_COMPOSITE_KEY_HPP
#define BOOST_MULTI_INDEX_COMPOSITE_KEY_HPP

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/multi_index/detail/access_specifier.hpp>
#include <boost/multi_index/detail/prevent_eti.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/aux_/nttp_decl.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp> 
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_same.hpp>
#include <functional>

/* A composite key stores n key extractors and "computes" the
 * result on a given value as a packed reference to the value and
 * the composite key itself. Actual invocations to the component
 * key extractors are lazily performed on comparison time.
 * As the other key extractors in Boost.MultiIndex, composite_key<T,...>
 * is  overloaded to work on chained pointers to T and reference_wrappers
 * of T.
 * composite_key_compare is overloaded to enable comparisons between
 * composite_key_results and tuples of values. Comparison is done
 * lexicographically on on the maximum common number of elements of the
 * operands. This allows searching for incomplete keys
 */

/* This user_definable macro limits the number of elements of a composite
 * key; useful for shortening resulting symbol names (MSVC++ 6.0, for
 * instance has problems coping with very long symbol names.)
 * NB: This cannot exceed the maximum number of arguments of
 * boost::tuple. In Boost 1.31, the limit is 10.
 */

#if !defined(BOOST_MULTI_INDEX_LIMIT_COMPOSITE_KEY_SIZE)
#if defined(BOOST_MSVC)&&(BOOST_MSVC<1300)
#define BOOST_MULTI_INDEX_LIMIT_COMPOSITE_KEY_SIZE 5
#else
#define BOOST_MULTI_INDEX_LIMIT_COMPOSITE_KEY_SIZE 10
#endif
#endif

/* maximum number of key extractors in a composite key */

#if BOOST_MULTI_INDEX_LIMIT_COMPOSITE_KEY_SIZE<10 /* max length of a tuple */
#define BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE \
  BOOST_MULTI_INDEX_LIMIT_COMPOSITE_KEY_SIZE
#else
#define BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE 10
#endif

/* BOOST_PP_ENUM of BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE elements */

#define BOOST_MULTI_INDEX_CK_ENUM(macro,data)                                \
  BOOST_PP_ENUM(BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE,macro,data)

/* BOOST_PP_ENUM_PARAMS of BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE elements */

#define BOOST_MULTI_INDEX_CK_ENUM_PARAMS(param)                              \
  BOOST_PP_ENUM_PARAMS(BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE,param)

/* if n==0 ->   text0
 * otherwise -> textn=tuples::null_type
 */

#define BOOST_MULTI_INDEX_CK_TEMPLATE_PARM(z,n,text)                         \
  typename BOOST_PP_CAT(text,n) BOOST_PP_EXPR_IF(n,=tuples::null_type)

/* const textn& kn=textn() */

#define BOOST_MULTI_INDEX_CK_CTOR_ARG(z,n,text)                              \
  const BOOST_PP_CAT(text,n)& BOOST_PP_CAT(k,n) = BOOST_PP_CAT(text,n)()

/* typename list(0)<list(1),n>::type */

#define BOOST_MULTI_INDEX_CK_APPLY_METAFUNCTION_N(z,n,list)                  \
  BOOST_DEDUCED_TYPENAME BOOST_PP_LIST_AT(list,0)<                           \
    BOOST_PP_LIST_AT(list,1),n                                               \
  >::type

namespace boost{

template<class T> class reference_wrapper; /* fwd decl. */

namespace multi_index{

namespace detail{

/* nth_composite_key_less<CompositeKey,N>:: yields std::less<result_type>
 * where result_type is the result_type of the nth key extractor of
 * CompositeKey. If N >= the length of CompositeKey, it yields
 * tuples::null_type.
 * Similar thing for nth_composite_key_greater.
 */

template<typename CompositeKey,BOOST_MPL_AUX_NTTP_DECL(int, N)>
struct nth_key_from_value
{
  typedef typename CompositeKey::key_extractor_tuple key_extractor_tuple;
  typedef typename prevent_eti<
    tuples::element<N,key_extractor_tuple>,
    typename mpl::eval_if_c<
      N<tuples::length<key_extractor_tuple>::value,
      tuples::element<N,key_extractor_tuple>,
      mpl::identity<tuples::null_type>
    >::type
  >::type                                            type;
};

template<typename KeyFromValue>
struct key_std_less
{
  typedef std::less<typename KeyFromValue::result_type> type;
};

template<>
struct key_std_less<tuples::null_type>
{
  typedef tuples::null_type type;
};

template<typename CompositeKey,BOOST_MPL_AUX_NTTP_DECL(int, N)>
struct nth_composite_key_less
{
  typedef typename nth_key_from_value<CompositeKey,N>::type key_from_value;
  typedef typename key_std_less<key_from_value>::type       type;
};

template<typename KeyFromValue>
struct key_std_greater
{
  typedef std::greater<typename KeyFromValue::result_type> type;
};

template<>
struct key_std_greater<tuples::null_type>
{
  typedef tuples::null_type type;
};

template<typename CompositeKey,BOOST_MPL_AUX_NTTP_DECL(int, N)>
struct nth_composite_key_greater
{
  typedef typename nth_key_from_value<CompositeKey,N>::type key_from_value;
  typedef typename key_std_greater<key_from_value>::type    type;
};

/* Metaprogramming machinery to compare composite_key_results between
 * them and with tuples of values.
 * equals_* computes equality of two tuple objects x,y with the same
 * length, defined as
 *
 *   xi==yi for all i in [0,...,min(length(x),length(y))).
 *
 * less_* accepts operands of different lenghts and computes the
 * following less-than relation:
 *
 *   !(xi<yi) && !(yi<xi) && xj<yj 
 *   for all i in [0,j) and some j in [0,...,min(length(x),length(y)).
 *
 * compare_* computes the same algorithm than less_*, but taking a tuple
 * of comparison predicates instead of operator<.
 */

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct equals_ckey_ckey; /* fwd decl. */

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct equals_ckey_ckey_terminal
{
  static bool compare(
    const KeyCons1&,const Value1&,
    const KeyCons2&,const Value2&)
  {
    return true;
  }
};

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct equals_ckey_ckey_normal
{
  static bool compare(
    const KeyCons1& c0,const Value1& v0,
    const KeyCons2& c1,const Value2& v1)
  {
    if(!(c0.get_head()(v0)==c1.get_head()(v1)))return false;
    return equals_ckey_ckey<
      BOOST_DEDUCED_TYPENAME KeyCons1::tail_type,Value1,
      BOOST_DEDUCED_TYPENAME KeyCons2::tail_type,Value2
    >::compare(c0.get_tail(),v0,c1.get_tail(),v1);
  }
};

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct equals_ckey_ckey:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons1,tuples::null_type>,
      is_same<KeyCons2,tuples::null_type>
    >,
    equals_ckey_ckey_terminal<KeyCons1,Value1,KeyCons2,Value2>,
    equals_ckey_ckey_normal<KeyCons1,Value1,KeyCons2,Value2>
  >::type
{
};

template<typename KeyCons,typename Value,typename ValCons>
struct equals_ckey_cval; /* fwd decl. */

template<typename KeyCons,typename Value,typename ValCons>
struct equals_ckey_cval_terminal
{
  static bool compare(const KeyCons&,const Value&,const ValCons&)
  {
    return true;
  }

  static bool compare(const ValCons&,const KeyCons&,const Value&)
  {
    return true;
  }
};

template<typename KeyCons,typename Value,typename ValCons>
struct equals_ckey_cval_normal
{
  static bool compare(const KeyCons& c,const Value& v,const ValCons& vc)
  {
    if(!(c.get_head()(v)==vc.get_head()))return false;
    return equals_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type
    >::compare(c.get_tail(),v,vc.get_tail());
  }

  static bool compare(const ValCons& vc,const KeyCons& c,const Value& v)
  {
    if(!(vc.get_head()==c.get_head()(v)))return false;
    return equals_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type
    >::compare(vc.get_tail(),c.get_tail(),v);
  }
};

template<typename KeyCons,typename Value,typename ValCons>
struct equals_ckey_cval:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons,tuples::null_type>,
      is_same<ValCons,tuples::null_type>
    >,
    equals_ckey_cval_terminal<KeyCons,Value,ValCons>,
    equals_ckey_cval_normal<KeyCons,Value, ValCons>
  >::type
{
};

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct less_ckey_ckey; /* fwd decl. */

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct less_ckey_ckey_terminal
{
  static bool compare(
    const KeyCons1&,const Value1&,const KeyCons2&,const Value2&)
  {
    return false;
  }
};

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct less_ckey_ckey_normal
{
  static bool compare(
    const KeyCons1& c0,const Value1& v0,
    const KeyCons2& c1,const Value2& v1)
  {
    if(c0.get_head()(v0)<c1.get_head()(v1))return true;
    if(c1.get_head()(v1)<c0.get_head()(v0))return false;
    return less_ckey_ckey<
      BOOST_DEDUCED_TYPENAME KeyCons1::tail_type,Value1,
      BOOST_DEDUCED_TYPENAME KeyCons2::tail_type,Value2
    >::compare(c0.get_tail(),v0,c1.get_tail(),v1);
  }
};

template<typename KeyCons1,typename Value1,typename KeyCons2,typename Value2>
struct less_ckey_ckey:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons1,tuples::null_type>,
      is_same<KeyCons2,tuples::null_type>
    >,
    less_ckey_ckey_terminal<KeyCons1,Value1,KeyCons2,Value2>,
    less_ckey_ckey_normal<KeyCons1,Value1,KeyCons2,Value2>
  >::type
{
};

template<typename KeyCons,typename Value,typename ValCons>
struct less_ckey_cval; /* fwd decl. */

template<typename KeyCons,typename Value,typename ValCons>
struct less_ckey_cval_terminal
{
  static bool compare(const KeyCons&,const Value&,const ValCons&)
  {
    return false;
  }

  static bool compare(const ValCons&,const KeyCons&,const Value&)
  {
    return false;
  }
};

template<typename KeyCons,typename Value,typename ValCons>
struct less_ckey_cval_normal
{
  static bool compare(const KeyCons& c,const Value& v,const ValCons& vc)
  {
    if(c.get_head()(v)<vc.get_head())return true;
    if(vc.get_head()<c.get_head()(v))return false;
    return less_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type
    >::compare(c.get_tail(),v,vc.get_tail());
  }

  static bool compare(const ValCons& vc,const KeyCons& c,const Value& v)
  {
    if(vc.get_head()<c.get_head()(v))return true;
    if(c.get_head()(v)<vc.get_head())return false;
    return less_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type
    >::compare(vc.get_tail(),c.get_tail(),v);
  }
};

template<typename KeyCons,typename Value,typename ValCons>
struct less_ckey_cval:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons,tuples::null_type>,
      is_same<ValCons,tuples::null_type>
    >,
    less_ckey_cval_terminal<KeyCons,Value,ValCons>,
    less_ckey_cval_normal<KeyCons,Value,ValCons>
  >::type
{
};

template
<
  typename KeyCons1,typename Value1,
  typename KeyCons2, typename Value2,
  typename CompareCons
>
struct compare_ckey_ckey; /* fwd decl. */

template
<
  typename KeyCons1,typename Value1,
  typename KeyCons2, typename Value2,
  typename CompareCons
>
struct compare_ckey_ckey_terminal
{
  static bool compare(
    const KeyCons1&,const Value1&,
    const KeyCons2&,const Value2&,
    const CompareCons&)
  {
    return false;
  }
};

template
<
  typename KeyCons1,typename Value1,
  typename KeyCons2, typename Value2,
  typename CompareCons
>
struct compare_ckey_ckey_normal
{
  static bool compare(
    const KeyCons1& c0,const Value1& v0,
    const KeyCons2& c1,const Value2& v1,
    const CompareCons& comp)
  {
    if(comp.get_head()(c0.get_head()(v0),c1.get_head()(v1)))return true;
    if(comp.get_head()(c1.get_head()(v1),c0.get_head()(v0)))return false;
    return compare_ckey_ckey<
      BOOST_DEDUCED_TYPENAME KeyCons1::tail_type,Value1,
      BOOST_DEDUCED_TYPENAME KeyCons2::tail_type,Value2,
      BOOST_DEDUCED_TYPENAME CompareCons::tail_type
    >::compare(c0.get_tail(),v0,c1.get_tail(),v1,comp.get_tail());
  }
};

template
<
  typename KeyCons1,typename Value1,
  typename KeyCons2, typename Value2,
  typename CompareCons
>
struct compare_ckey_ckey:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons1,tuples::null_type>,
      is_same<KeyCons2,tuples::null_type>
    >,
    compare_ckey_ckey_terminal<KeyCons1,Value1,KeyCons2,Value2,CompareCons>,
    compare_ckey_ckey_normal<KeyCons1,Value1,KeyCons2,Value2,CompareCons>
  >::type
{
};

template
<
  typename KeyCons,typename Value,
  typename ValCons,typename CompareCons
>
struct compare_ckey_cval; /* fwd decl. */

template
<
  typename KeyCons,typename Value,
  typename ValCons,typename CompareCons
>
struct compare_ckey_cval_terminal
{
  static bool compare(
    const KeyCons&,const Value&,const ValCons&,const CompareCons&)
  {
    return false;
  }

  static bool compare(
    const ValCons&,const KeyCons&,const Value&,const CompareCons&)
  {
    return false;
  }
};

template
<
  typename KeyCons,typename Value,
  typename ValCons,typename CompareCons
>
struct compare_ckey_cval_normal
{
  static bool compare(
    const KeyCons& c,const Value& v,const ValCons& vc,
    const CompareCons& comp)
  {
    if(comp.get_head()(c.get_head()(v),vc.get_head()))return true;
    if(comp.get_head()(vc.get_head(),c.get_head()(v)))return false;
    return compare_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type,
      BOOST_DEDUCED_TYPENAME CompareCons::tail_type
    >::compare(c.get_tail(),v,vc.get_tail(),comp.get_tail());
  }

  static bool compare(
    const ValCons& vc,const KeyCons& c,const Value& v,
    const CompareCons& comp)
  {
    if(comp.get_head()(vc.get_head(),c.get_head()(v)))return true;
    if(comp.get_head()(c.get_head()(v),vc.get_head()))return false;
    return compare_ckey_cval<
      BOOST_DEDUCED_TYPENAME KeyCons::tail_type,Value,
      BOOST_DEDUCED_TYPENAME ValCons::tail_type,
      BOOST_DEDUCED_TYPENAME CompareCons::tail_type
    >::compare(vc.get_tail(),c.get_tail(),v,comp.get_tail());
  }
};

template
<
  typename KeyCons,typename Value,
  typename ValCons,typename CompareCons
>
struct compare_ckey_cval:
  mpl::if_<
    mpl::or_<
      is_same<KeyCons,tuples::null_type>,
      is_same<ValCons,tuples::null_type>
    >,
    compare_ckey_cval_terminal<KeyCons,Value,ValCons,CompareCons>,
    compare_ckey_cval_normal<KeyCons,Value,ValCons,CompareCons>
  >::type
{
};

} /* namespace multi_index::detail */

/* composite_key_result */

template<typename CompositeKey>
struct composite_key_result
{
  typedef CompositeKey                            composite_key_type;
  typedef typename composite_key_type::value_type value_type;

  composite_key_result(
    const composite_key_type& composite_key_,const value_type& value_):
    composite_key(composite_key_),value(value_)
  {}

  const composite_key_type& composite_key;
  const value_type&         value;
};

/* composite_key */

/* NB. Some overloads of operator() have an extra dummy parameter int=0.
 * This disambiguator serves several purposes:
 *  - Without it, MSVC++ 6.0 incorrectly regards some overloads as
 *    specializations of a previous member function template.
 *  - MSVC++ 6.0/7.0 seem to incorrectly treat some different memfuns
 *    as if they have the same signature.
 *  - If remove_const is broken due to lack of PTS, int=0 avoids the
 *    declaration of memfuns with identical signature.
 */

template<
  typename Value,
  BOOST_MULTI_INDEX_CK_ENUM(BOOST_MULTI_INDEX_CK_TEMPLATE_PARM,KeyFromValue)
>
struct composite_key:
  private tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(KeyFromValue)>
{
private:
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(KeyFromValue)> super;

public:
  typedef super                               key_extractor_tuple;
  typedef Value                               value_type;
  typedef composite_key_result<composite_key> result_type;

  composite_key(
    BOOST_MULTI_INDEX_CK_ENUM(BOOST_MULTI_INDEX_CK_CTOR_ARG,KeyFromValue)):
    super(BOOST_MULTI_INDEX_CK_ENUM_PARAMS(k))
  {}

  composite_key(const key_extractor_tuple& x):super(x){}

  const key_extractor_tuple& key_extractors()const{return *this;}
  key_extractor_tuple&       key_extractors(){return *this;}

  template<typename ChainedPtr>
  result_type operator()(const ChainedPtr& x)const
  {
    return operator()(*x);
  }

  result_type operator()(const value_type& x)const
  {
    return result_type(*this,x);
  }

  result_type operator()(const reference_wrapper<const value_type>& x)const
  {
    return result_type(*this,x.get());
  }

  result_type operator()(const reference_wrapper<value_type>& x,int=0)const
  {
    return result_type(*this,x.get());
  }
};

/* comparison operators */

/* == */

template<typename CompositeKey1,typename CompositeKey2>
inline bool operator==(
  const composite_key_result<CompositeKey1>& x,
  const composite_key_result<CompositeKey2>& y)
{
  typedef typename CompositeKey1::key_extractor_tuple key_extractor_tuple1;
  typedef typename CompositeKey1::value_type          value_type1;
  typedef typename CompositeKey2::key_extractor_tuple key_extractor_tuple2;
  typedef typename CompositeKey2::value_type          value_type2;

  BOOST_STATIC_ASSERT(
    tuples::length<key_extractor_tuple1>::value==
    tuples::length<key_extractor_tuple2>::value);

  return detail::equals_ckey_ckey<
    key_extractor_tuple1,value_type1,
    key_extractor_tuple2,value_type2
  >::compare(
    x.composite_key.key_extractors(),x.value,
    y.composite_key.key_extractors(),y.value);
}

template<
  typename CompositeKey,
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value)
>
inline bool operator==(
  const composite_key_result<CompositeKey>& x,
  const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& y)
{
  typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
  typedef typename CompositeKey::value_type              value_type;
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;
  
  BOOST_STATIC_ASSERT(
    tuples::length<key_extractor_tuple>::value==
    tuples::length<key_tuple>::value);

  return detail::equals_ckey_cval<key_extractor_tuple,value_type,key_tuple>::
    compare(x.composite_key.key_extractors(),x.value,y);
}

template
<
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value),
  typename CompositeKey
>
inline bool operator==(
  const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& x,
  const composite_key_result<CompositeKey>& y)
{
  typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
  typedef typename CompositeKey::value_type              value_type;
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;
  
  BOOST_STATIC_ASSERT(
    tuples::length<key_extractor_tuple>::value==
    tuples::length<key_tuple>::value);

  return detail::equals_ckey_cval<key_extractor_tuple,value_type,key_tuple>::
    compare(x,y.composite_key.key_extractors(),y.value);
}

/* < */

template<typename CompositeKey1,typename CompositeKey2>
inline bool operator<(
  const composite_key_result<CompositeKey1>& x,
  const composite_key_result<CompositeKey2>& y)
{
  typedef typename CompositeKey1::key_extractor_tuple key_extractor_tuple1;
  typedef typename CompositeKey1::value_type          value_type1;
  typedef typename CompositeKey2::key_extractor_tuple key_extractor_tuple2;
  typedef typename CompositeKey2::value_type          value_type2;

  return detail::less_ckey_ckey<
   key_extractor_tuple1,value_type1,
   key_extractor_tuple2,value_type2
  >::compare(
    x.composite_key.key_extractors(),x.value,
    y.composite_key.key_extractors(),y.value);
}

template
<
  typename CompositeKey,
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value)
>
inline bool operator<(
  const composite_key_result<CompositeKey>& x,
  const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& y)
{
  typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
  typedef typename CompositeKey::value_type              value_type;
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;
  
  return detail::less_ckey_cval<key_extractor_tuple,value_type,key_tuple>::
    compare(x.composite_key.key_extractors(),x.value,y);
}

template
<
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value),
  typename CompositeKey
>
inline bool operator<(
  const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& x,
  const composite_key_result<CompositeKey>& y)
{
  typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
  typedef typename CompositeKey::value_type              value_type;
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;
  
  return detail::less_ckey_cval<key_extractor_tuple,value_type,key_tuple>::
    compare(x,y.composite_key.key_extractors(),y.value);
}

/* rest of comparison operators */

#define BOOST_MULTI_INDEX_CK_COMPLETE_COMP_OPS(t1,t2,a1,a2)                  \
template<t1,t2> inline bool operator!=(const a1& x,const a2& y)              \
{                                                                            \
  return !(x==y);                                                            \
}                                                                            \
                                                                             \
template<t1,t2> inline bool operator>(const a1& x,const a2& y)               \
{                                                                            \
  return y<x;                                                                \
}                                                                            \
                                                                             \
template<t1,t2> inline bool operator>=(const a1& x,const a2& y)              \
{                                                                            \
  return !(x<y);                                                             \
}                                                                            \
                                                                             \
template<t1,t2> inline bool operator<=(const a1& x,const a2& y)              \
{                                                                            \
  return !(y<x);                                                             \
}

BOOST_MULTI_INDEX_CK_COMPLETE_COMP_OPS(
  typename CompositeKey1,
  typename CompositeKey2,
  composite_key_result<CompositeKey1>,
  composite_key_result<CompositeKey2>
)

BOOST_MULTI_INDEX_CK_COMPLETE_COMP_OPS(
  typename CompositeKey,
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value),
  composite_key_result<CompositeKey>,
  tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>
)

BOOST_MULTI_INDEX_CK_COMPLETE_COMP_OPS(
  BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value),
  typename CompositeKey,
  tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>,
  composite_key_result<CompositeKey>
)

/* composite_key_compare */

template
<
  BOOST_MULTI_INDEX_CK_ENUM(BOOST_MULTI_INDEX_CK_TEMPLATE_PARM,Compare)
>
struct composite_key_compare:
  private tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Compare)>
{
private:
  typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Compare)> super;

public:
  typedef super key_comp_tuple;

  composite_key_compare(
    BOOST_MULTI_INDEX_CK_ENUM(BOOST_MULTI_INDEX_CK_CTOR_ARG,Compare)):
    super(BOOST_MULTI_INDEX_CK_ENUM_PARAMS(k))
  {}

  composite_key_compare(const key_comp_tuple& x):super(x){}

  const key_comp_tuple& key_comps()const{return *this;}
  key_comp_tuple&       key_comps(){return *this;}

  template<typename CompositeKey1,typename CompositeKey2>
  bool operator()(
    const composite_key_result<CompositeKey1> & x,
    const composite_key_result<CompositeKey2> & y)const
  {
    typedef typename CompositeKey1::key_extractor_tuple key_extractor_tuple1;
    typedef typename CompositeKey1::value_type          value_type1;
    typedef typename CompositeKey2::key_extractor_tuple key_extractor_tuple2;
    typedef typename CompositeKey2::value_type          value_type2;

    BOOST_STATIC_ASSERT(
      tuples::length<key_extractor_tuple1>::value<=
      tuples::length<key_comp_tuple>::value||
      tuples::length<key_extractor_tuple2>::value<=
      tuples::length<key_comp_tuple>::value);

    return detail::compare_ckey_ckey<
      key_extractor_tuple1,value_type1,
      key_extractor_tuple2,value_type2,
      key_comp_tuple
    >::compare(
      x.composite_key.key_extractors(),x.value,
      y.composite_key.key_extractors(),y.value,
      key_comps());
  }
  
  template
  <
    typename CompositeKey,
    BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value)
  >
  bool operator()(
    const composite_key_result<CompositeKey>& x,
    const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& y)const
  {
    typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
    typedef typename CompositeKey::value_type              value_type;
    typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;

    BOOST_STATIC_ASSERT(
      tuples::length<key_extractor_tuple>::value<=
      tuples::length<key_comp_tuple>::value||
      tuples::length<key_tuple>::value<=
      tuples::length<key_comp_tuple>::value);

    return detail::compare_ckey_cval<
      key_extractor_tuple,value_type,
      key_tuple,key_comp_tuple
    >::compare(x.composite_key.key_extractors(),x.value,y,key_comps());
  }

  template
  <
    BOOST_MULTI_INDEX_CK_ENUM_PARAMS(typename Value),
    typename CompositeKey
  >
  bool operator()(
    const tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)>& x,
    const composite_key_result<CompositeKey>& y)const
  {
    typedef typename CompositeKey::key_extractor_tuple     key_extractor_tuple;
    typedef typename CompositeKey::value_type              value_type;
    typedef tuple<BOOST_MULTI_INDEX_CK_ENUM_PARAMS(Value)> key_tuple;

    BOOST_STATIC_ASSERT(
      tuples::length<key_tuple>::value<=
      tuples::length<key_comp_tuple>::value||
      tuples::length<key_extractor_tuple>::value<=
      tuples::length<key_comp_tuple>::value);

    return detail::compare_ckey_cval<
      key_extractor_tuple,value_type,
      key_tuple,key_comp_tuple
    >::compare(x,y.composite_key.key_extractors(),y.value,key_comps());
  }
};

/* composite_key_compare_less is merely a composite_key_compare
 * instantiation with the corresponding std::less<> comparison
 * predicates for each key extractor. Useful as a substitute for
 * std::less<CompositeKey::result_type> when the compiler does not
 * support partial specialization.
 * Same with composite_key_compare_greater.
 */

#define BOOST_MULTI_INDEX_CK_RESULT_LESS_SUPER                               \
composite_key_compare<                                                       \
    BOOST_MULTI_INDEX_CK_ENUM(                                               \
      BOOST_MULTI_INDEX_CK_APPLY_METAFUNCTION_N,                             \
      /* the argument is a PP list */                                        \
      (detail::nth_composite_key_less,                                       \
        (BOOST_DEDUCED_TYPENAME CompositeKeyResult::composite_key_type,      \
          BOOST_PP_NIL)))                                                    \
  >

template<typename CompositeKeyResult>
struct composite_key_result_less:
BOOST_MULTI_INDEX_PRIVATE_IF_USING_DECL_FOR_TEMPL_FUNCTIONS
BOOST_MULTI_INDEX_CK_RESULT_LESS_SUPER
{
private:
  typedef BOOST_MULTI_INDEX_CK_RESULT_LESS_SUPER super;

public:
  typedef CompositeKeyResult  first_argument_type;
  typedef first_argument_type second_argument_type;
  typedef bool                result_type;

  using super::operator();
};

#define BOOST_MULTI_INDEX_CK_RESULT_GREATER_SUPER                            \
composite_key_compare<                                                       \
    BOOST_MULTI_INDEX_CK_ENUM(                                               \
      BOOST_MULTI_INDEX_CK_APPLY_METAFUNCTION_N,                             \
      /* the argument is a PP list */                                        \
      (detail::nth_composite_key_greater,                                    \
        (BOOST_DEDUCED_TYPENAME CompositeKeyResult::composite_key_type,      \
          BOOST_PP_NIL)))                                                    \
  >

template<typename CompositeKeyResult>
struct composite_key_result_greater:
BOOST_MULTI_INDEX_PRIVATE_IF_USING_DECL_FOR_TEMPL_FUNCTIONS
BOOST_MULTI_INDEX_CK_RESULT_GREATER_SUPER
{
private:
  typedef BOOST_MULTI_INDEX_CK_RESULT_GREATER_SUPER super;

public:
  typedef CompositeKeyResult  first_argument_type;
  typedef first_argument_type second_argument_type;
  typedef bool                result_type;

  using super::operator();
};

} /* namespace multi_index */

} /* namespace boost */

/* Specialization of std::less and std::greater for composite_key_results
 * enabling comparison with tuples of values.
 */

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
namespace std{

template<typename CompositeKey>
struct less<boost::multi_index::composite_key_result<CompositeKey> >:
  boost::multi_index::composite_key_result_less<
    boost::multi_index::composite_key_result<CompositeKey>
  >
{
};

template<typename CompositeKey>
struct greater<boost::multi_index::composite_key_result<CompositeKey> >:
  boost::multi_index::composite_key_result_greater<
    boost::multi_index::composite_key_result<CompositeKey>
  >
{
};

} /* namespace std */
#endif

#undef BOOST_MULTI_INDEX_CK_RESULT_LESS_SUPER
#undef BOOST_MULTI_INDEX_CK_RESULT_GREATER_SUPER
#undef BOOST_MULTI_INDEX_CK_COMPLETE_COMP_OPS
#undef BOOST_MULTI_INDEX_CK_APPLY_METAFUNCTION_N
#undef BOOST_MULTI_INDEX_CK_CTOR_ARG
#undef BOOST_MULTI_INDEX_CK_TEMPLATE_PARM
#undef BOOST_MULTI_INDEX_CK_ENUM_PARAMS
#undef BOOST_MULTI_INDEX_CK_ENUM
#undef BOOST_MULTI_INDEX_COMPOSITE_KEY_SIZE

#endif
