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
# ifndef BOOST_PREPROCESSOR_LIST_TO_TUPLE_HPP
# define BOOST_PREPROCESSOR_LIST_TO_TUPLE_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/list/enum.hpp>
#
# /* BOOST_PP_LIST_TO_TUPLE */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_LIST_TO_TUPLE(list) (BOOST_PP_LIST_ENUM(list))
# else
#    define BOOST_PP_LIST_TO_TUPLE(list) BOOST_PP_LIST_TO_TUPLE_I(list)
#    define BOOST_PP_LIST_TO_TUPLE_I(list) (BOOST_PP_LIST_ENUM(list))
# endif
#
# /* BOOST_PP_LIST_TO_TUPLE_R */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_LIST_TO_TUPLE_R(r, list) (BOOST_PP_LIST_ENUM_R(r, list))
# else
#    define BOOST_PP_LIST_TO_TUPLE_R(r, list) BOOST_PP_LIST_TO_TUPLE_R_I(r, list)
#    define BOOST_PP_LIST_TO_TUPLE_R_I(r, list) (BOOST_PP_LIST_ENUM_R(r, list))
# endif
#
# endif
