
#ifndef BOOST_MPL_AUX_PREPROCESSOR_IS_SEQ_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_IS_SEQ_HPP_INCLUDED

// + file: boost/mpl/aux_/preprocessor/is_seq.hpp
// + last modified: 03/may/03

// Copyright (c) 2003
// Paul Mensonides, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/preprocessor/seq/size.hpp"
#include "boost/preprocessor/arithmetic/dec.hpp"
#include "boost/preprocessor/punctuation/paren.hpp"
#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/config/config.hpp"

// returns 1 if 'seq' is a PP-sequence, 0 otherwise:
//
//   BOOST_PP_ASSERT( BOOST_PP_NOT( BOOST_MPL_PP_IS_SEQ( int ) ) )
//   BOOST_PP_ASSERT( BOOST_MPL_PP_IS_SEQ( (int) ) )
//   BOOST_PP_ASSERT( BOOST_MPL_PP_IS_SEQ( (1)(2) ) )

#if BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_BCC()

#   define BOOST_MPL_PP_IS_SEQ(seq) BOOST_PP_DEC( BOOST_PP_SEQ_SIZE( BOOST_MPL_PP_IS_SEQ_(seq) ) )
#   define BOOST_MPL_PP_IS_SEQ_(seq) BOOST_MPL_PP_IS_SEQ_SEQ_( BOOST_MPL_PP_IS_SEQ_SPLIT_ seq )
#   define BOOST_MPL_PP_IS_SEQ_SEQ_(x) (x)
#   define BOOST_MPL_PP_IS_SEQ_SPLIT_(unused) unused)((unused)

#else

#   if BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#       define BOOST_MPL_PP_IS_SEQ(seq) BOOST_MPL_PP_IS_SEQ_MWCC_((seq))
#       define BOOST_MPL_PP_IS_SEQ_MWCC_(args) BOOST_MPL_PP_IS_SEQ_ ## args
#   else
#       define BOOST_MPL_PP_IS_SEQ(seq) BOOST_MPL_PP_IS_SEQ_(seq)
#   endif

#   define BOOST_MPL_PP_IS_SEQ_(seq) BOOST_PP_CAT(BOOST_MPL_PP_IS_SEQ_, BOOST_MPL_PP_IS_SEQ_0 seq BOOST_PP_RPAREN())
#   define BOOST_MPL_PP_IS_SEQ_0(x) BOOST_MPL_PP_IS_SEQ_1(x
#   define BOOST_MPL_PP_IS_SEQ_ALWAYS_0(unused) 0
#   define BOOST_MPL_PP_IS_SEQ_BOOST_MPL_PP_IS_SEQ_0 BOOST_MPL_PP_IS_SEQ_ALWAYS_0(
#   define BOOST_MPL_PP_IS_SEQ_BOOST_MPL_PP_IS_SEQ_1(unused) 1

#endif

#endif // BOOST_MPL_AUX_PREPROCESSOR_IS_SEQ_HPP_INCLUDED
