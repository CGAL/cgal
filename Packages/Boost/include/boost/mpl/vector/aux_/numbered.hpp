//-----------------------------------------------------------------------------
// boost mpl/vector/aux_/numbered.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

// no include guards, the header is intended for multiple inclusion!

#if defined(BOOST_PP_IS_ITERATING)

#include "boost/preprocessor/enum_params.hpp"
#include "boost/preprocessor/enum_shifted_params.hpp"
#include "boost/preprocessor/comma_if.hpp"
#include "boost/preprocessor/repeat.hpp"
#include "boost/preprocessor/dec.hpp"
#include "boost/preprocessor/cat.hpp"

#define i BOOST_PP_FRAME_ITERATION(1)

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)

#   define MPL_AUX_VECTOR_TAIL(vector, i, T) \
    BOOST_PP_CAT(vector,BOOST_PP_DEC(i))< \
      BOOST_PP_ENUM_SHIFTED_PARAMS(i, T) \
    > \
    /**/

#if i > 0
template<
      BOOST_PP_ENUM_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(vector,i)
    : vector_node<
          i
        , T0
        , MPL_AUX_VECTOR_TAIL(vector,i,T)
        >
{
};
#endif

#   undef MPL_AUX_VECTOR_TAIL

#else // "brute force" implementation

#   if i > 0

template<
      BOOST_PP_ENUM_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(vector,i)
{
    typedef aux::vector_tag<i> tag;
    typedef BOOST_PP_CAT(vector,i) type;

#   define AUX_VECTOR_ITEM(unused, i, unused2) \
    typedef BOOST_PP_CAT(T,i) BOOST_PP_CAT(item,i); \
    /**/

    BOOST_PP_REPEAT_1(i, AUX_VECTOR_ITEM, unused)
#   undef AUX_VECTOR_ITEM
    typedef void_ BOOST_PP_CAT(item,i);
    typedef BOOST_PP_CAT(T,BOOST_PP_DEC(i)) back;

    // Borland forces us to use |type| here (instead of the class name)
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,i> > end;
};

template<>
struct push_front_traits< aux::vector_tag<BOOST_PP_DEC(i)> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef BOOST_PP_CAT(vector,i)<
              T
              BOOST_PP_COMMA_IF(BOOST_PP_DEC(i))
              BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(i), typename Vector::item)
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
    {
        typedef BOOST_PP_CAT(vector,BOOST_PP_DEC(i))<
              BOOST_PP_ENUM_SHIFTED_PARAMS(i, typename Vector::item)
            > type;
    };
};

#   endif // i > 0

#   if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
    && !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)

template< typename V >
struct vector_item<V,i>
{
    typedef typename V::BOOST_PP_CAT(item,i) type;
};

#   else

namespace aux {
template<> struct vector_item_impl<i>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::BOOST_PP_CAT(item,i) type;
    };
};
}

template<>
struct at_traits< aux::vector_tag<i> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

#if i > 0
template<>
struct front_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};
#endif

template<>
struct size_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
        : integral_c<long,i>
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag<i> >
    : size_traits< aux::vector_tag<i> >
{
};

template<>
struct clear_traits< aux::vector_tag<i> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#endif // BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL

#undef i

#endif // BOOST_PP_IS_ITERATING
