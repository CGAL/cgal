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
# ifndef BOOST_PREPROCESSOR_LIST_ENUM_HPP
# define BOOST_PREPROCESSOR_LIST_ENUM_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/list/for_each_i.hpp>
# include <boost/preprocessor/punctuation/comma_if.hpp>
#
# /* BOOST_PP_LIST_ENUM */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_LIST_ENUM(list) BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_LIST_ENUM_O, BOOST_PP_NIL, list)
# else
#    define BOOST_PP_LIST_ENUM(list) BOOST_PP_LIST_ENUM_I(list)
#    define BOOST_PP_LIST_ENUM_I(list) BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_LIST_ENUM_O, BOOST_PP_NIL, list)
# endif
#
# define BOOST_PP_LIST_ENUM_O(r, _, i, elem) BOOST_PP_COMMA_IF(i) elem
#
# /* BOOST_PP_LIST_ENUM_R */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_LIST_ENUM_R(r, list) BOOST_PP_LIST_FOR_EACH_I_R(r, BOOST_PP_LIST_ENUM_O, BOOST_PP_NIL, list)
# else
#    define BOOST_PP_LIST_ENUM_R(r, list) BOOST_PP_LIST_ENUM_R_I(r, list)
#    define BOOST_PP_LIST_ENUM_R_I(r, list) BOOST_PP_LIST_FOR_EACH_I_R(r, BOOST_PP_LIST_ENUM_O, BOOST_PP_NIL, list)
# endif
#
# endif
