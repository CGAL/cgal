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
# ifndef BOOST_PREPROCESSOR_SEQ_REVERSE_HPP
# define BOOST_PREPROCESSOR_SEQ_REVERSE_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/seq/fold_left.hpp>
#
# /* BOOST_PP_SEQ_REVERSE */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_REVERSE(seq) BOOST_PP_SEQ_FOLD_LEFT(BOOST_PP_SEQ_REVERSE_O, BOOST_PP_EMPTY, seq)()
# else
#    define BOOST_PP_SEQ_REVERSE(seq) BOOST_PP_SEQ_REVERSE_I(seq)
#    define BOOST_PP_SEQ_REVERSE_I(seq) BOOST_PP_SEQ_FOLD_LEFT(BOOST_PP_SEQ_REVERSE_O, BOOST_PP_EMPTY, seq)()
# endif
#
# define BOOST_PP_SEQ_REVERSE_O(s, state, elem) (elem) state
#
# /* BOOST_PP_SEQ_REVERSE_S */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_REVERSE_S(s, seq) BOOST_PP_SEQ_FOLD_LEFT_ ## s(BOOST_PP_SEQ_REVERSE_O, BOOST_PP_EMPTY, seq)()
# else
#    define BOOST_PP_SEQ_REVERSE_S(s, seq) BOOST_PP_SEQ_REVERSE_S_I(s, seq)
#    define BOOST_PP_SEQ_REVERSE_S_I(s, seq) BOOST_PP_SEQ_FOLD_LEFT_ ## s(BOOST_PP_SEQ_REVERSE_O, BOOST_PP_EMPTY, seq)()
# endif
#
# endif
