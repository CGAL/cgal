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
# ifndef BOOST_PREPROCESSOR_ARRAY_ELEM_HPP
# define BOOST_PREPROCESSOR_ARRAY_ELEM_HPP
#
# include <boost/preprocessor/array/data.hpp>
# include <boost/preprocessor/array/size.hpp>
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/tuple/elem.hpp>
#
# /* BOOST_PP_ARRAY_ELEM */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_ARRAY_ELEM(i, array) BOOST_PP_TUPLE_ELEM(BOOST_PP_ARRAY_SIZE(array), i, BOOST_PP_ARRAY_DATA(array))
# else
#    define BOOST_PP_ARRAY_ELEM(i, array) BOOST_PP_ARRAY_ELEM_I(i, array)
#    define BOOST_PP_ARRAY_ELEM_I(i, array) BOOST_PP_TUPLE_ELEM(BOOST_PP_ARRAY_SIZE(array), i, BOOST_PP_ARRAY_DATA(array))
# endif
#
# endif
