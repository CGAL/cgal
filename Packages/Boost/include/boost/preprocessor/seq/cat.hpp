# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2002.  Permission to copy, use,        *
#  *     modify, sell, and distribute this software is granted provided       *
#  *     this copyright notice appears in all copies.  This software is       *
#  *     provided "as is" without express or implied warranty, and with       *
#  *     no claim at to its suitability for any purpose.                      *
#  *                                                                          *
#  ************************************************************************** */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_SEQ_CAT_HPP
# define BOOST_PREPROCESSOR_SEQ_CAT_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/seq/fold_left.hpp>
# include <boost/preprocessor/seq/seq.hpp>
#
# /* BOOST_PP_SEQ_CAT */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_CAT(seq) BOOST_PP_SEQ_FOLD_LEFT(BOOST_PP_SEQ_CAT_O, BOOST_PP_SEQ_HEAD(seq), BOOST_PP_SEQ_TAIL(seq))
# else
#    define BOOST_PP_SEQ_CAT(seq) BOOST_PP_SEQ_CAT_I(seq)
#    define BOOST_PP_SEQ_CAT_I(seq) BOOST_PP_SEQ_FOLD_LEFT(BOOST_PP_SEQ_CAT_O, BOOST_PP_SEQ_HEAD(seq), BOOST_PP_SEQ_TAIL(seq))
# endif
#
# define BOOST_PP_SEQ_CAT_O(s, st, elem) BOOST_PP_SEQ_CAT_O_I(st, elem)
# define BOOST_PP_SEQ_CAT_O_I(a, b) a ## b
#
# /* BOOST_PP_SEQ_CAT_S */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_CAT_S(s, seq) BOOST_PP_SEQ_FOLD_LEFT_ ## s(BOOST_PP_SEQ_CAT_O, BOOST_PP_SEQ_HEAD(seq), BOOST_PP_SEQ_TAIL(seq))
# else
#    define BOOST_PP_SEQ_CAT_S(s, seq) BOOST_PP_SEQ_CAT_S_I(s, seq)
#    define BOOST_PP_SEQ_CAT_S_I(s, seq) BOOST_PP_SEQ_FOLD_LEFT_ ## s(BOOST_PP_SEQ_CAT_O, BOOST_PP_SEQ_HEAD(seq), BOOST_PP_SEQ_TAIL(seq))
# endif
#
# endif
