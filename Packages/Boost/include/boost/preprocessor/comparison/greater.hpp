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
# ifndef BOOST_PREPROCESSOR_COMPARISON_GREATER_HPP
# define BOOST_PREPROCESSOR_COMPARISON_GREATER_HPP
#
# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_GREATER */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_GREATER(x, y) BOOST_PP_LESS(y, x)
# else
#    define BOOST_PP_GREATER(x, y) BOOST_PP_GREATER_I(x, y)
#    define BOOST_PP_GREATER_I(x, y) BOOST_PP_LESS(y, x)
# endif
#
# /* BOOST_PP_GREATER_D */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_GREATER_D(d, x, y) BOOST_PP_LESS_D(d, y, x)
# else
#    define BOOST_PP_GREATER_D(d, x, y) BOOST_PP_GREATER_D_I(d, x, y)
#    define BOOST_PP_GREATER_D_I(d, x, y) BOOST_PP_LESS_D(d, y, x)
# endif
#
# endif
