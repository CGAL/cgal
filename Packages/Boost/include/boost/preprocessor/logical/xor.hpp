# /* Copyright (C) 2001
#  * Housemarque Oy
#  * http://www.housemarque.com
#  *
#  * Permission to copy, use, modify, sell and distribute this software is
#  * granted provided this copyright notice appears in all copies. This
#  * software is provided "as is" without express or implied warranty, and
#  * with no claim as to its suitability for any purpose.
#  */
#
# /* Revised by Paul Mensonides (2002) */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_LOGICAL_XOR_HPP
# define BOOST_PREPROCESSOR_LOGICAL_XOR_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/logical/bool.hpp>
# include <boost/preprocessor/logical/bitxor.hpp>
#
# /* BOOST_PP_XOR */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_XOR(p, q) BOOST_PP_BITXOR(BOOST_PP_BOOL(p), BOOST_PP_BOOL(q))
# else
#    define BOOST_PP_XOR(p, q) BOOST_PP_XOR_I(p, q)
#    define BOOST_PP_XOR_I(p, q) BOOST_PP_BITXOR(BOOST_PP_BOOL(p), BOOST_PP_BOOL(q))
# endif
#
# endif
