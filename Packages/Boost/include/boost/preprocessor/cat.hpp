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
# ifndef BOOST_PREPROCESSOR_CAT_HPP
# define BOOST_PREPROCESSOR_CAT_HPP
#
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_CAT */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#    define BOOST_PP_CAT(a, b) BOOST_PP_CAT_I(a, b)
# else
#    define BOOST_PP_CAT(a, b) BOOST_PP_CAT_OO((a, b))
#    define BOOST_PP_CAT_OO(par) BOOST_PP_CAT_I ## par
# endif
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MSVC()
#    define BOOST_PP_CAT_I(a, b) a ## b
# else
#    define BOOST_PP_CAT_I(a, b) BOOST_PP_CAT_II(a ## b)
#    define BOOST_PP_CAT_II(res) res
# endif
#
# endif
