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
# ifndef BOOST_PREPROCESSOR_DETAIL_IS_NULLARY_HPP
# define BOOST_PREPROCESSOR_DETAIL_IS_NULLARY_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/detail/check.hpp>
#
# /* BOOST_PP_IS_NULLARY */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_IS_NULLARY(x) BOOST_PP_CHECK(x, BOOST_PP_IS_NULLARY_CHECK)
# else
#    define BOOST_PP_IS_NULLARY(x) BOOST_PP_IS_NULLARY_I(x)
#    define BOOST_PP_IS_NULLARY_I(x) BOOST_PP_CHECK(x, BOOST_PP_IS_NULLARY_CHECK)
# endif
#
# define BOOST_PP_IS_NULLARY_CHECK() 1
# define BOOST_PP_CHECK_RESULT_BOOST_PP_IS_NULLARY_CHECK 0, BOOST_PP_NIL
#
# endif
