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
# ifndef BOOST_PREPROCESSOR_CONTROL_EXPR_IF_HPP
# define BOOST_PREPROCESSOR_CONTROL_EXPR_IF_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/control/expr_iif.hpp>
# include <boost/preprocessor/logical/bool.hpp>
#
# /* BOOST_PP_EXPR_IF */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_EXPR_IF(cond, expr) BOOST_PP_EXPR_IIF(BOOST_PP_BOOL(cond), expr)
# else
#    define BOOST_PP_EXPR_IF(cond, expr) BOOST_PP_EXPR_IF_I(cond, expr)
#    define BOOST_PP_EXPR_IF_I(cond, expr) BOOST_PP_EXPR_IIF(BOOST_PP_BOOL(cond), expr)
# endif
#
# endif
