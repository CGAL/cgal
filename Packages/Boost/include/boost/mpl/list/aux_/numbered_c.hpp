//-----------------------------------------------------------------------------
// boost mpl/list/aux_/numbered_c.hpp header file
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
#include "boost/preprocessor/dec.hpp"
#include "boost/preprocessor/cat.hpp"

#define i BOOST_PP_FRAME_ITERATION(1)

#if i == 1

template<
      typename T
    , BOOST_PP_ENUM_PARAMS(i, T C)
    >
struct list1_c
    : list_node<
          integral_c<long,1>
        , integral_c<T,C0>
        , null_node
        >
{
    typedef list1_c type;
};

#else

#   define MPL_AUX_LIST_C_TAIL(list, i, C) \
    BOOST_PP_CAT(BOOST_PP_CAT(list,BOOST_PP_DEC(i)),_c)<T, \
      BOOST_PP_ENUM_SHIFTED_PARAMS(i, C) \
    > \
    /**/
    
template<
      typename T
    , BOOST_PP_ENUM_PARAMS(i, T C)
    >
struct BOOST_PP_CAT(BOOST_PP_CAT(list,i),_c)
    : list_node<
          integral_c<long,i>
        , integral_c<T,C0>
        , MPL_AUX_LIST_C_TAIL(list,i,C)
        >
{
    typedef BOOST_PP_CAT(BOOST_PP_CAT(list,i),_c) type;
};

#   undef MPL_AUX_LIST_C_TAIL

#endif // i == 1

#undef i

#endif // BOOST_PP_IS_ITERATING
