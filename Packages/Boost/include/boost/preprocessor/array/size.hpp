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
# ifndef BOOST_PREPROCESSOR_ARRAY_SIZE_HPP
# define BOOST_PREPROCESSOR_ARRAY_SIZE_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/tuple/elem.hpp>
#
# /* BOOST_PP_ARRAY_SIZE */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_ARRAY_SIZE(array) BOOST_PP_TUPLE_ELEM(2, 0, array)
# else
#    define BOOST_PP_ARRAY_SIZE(array) BOOST_PP_ARRAY_SIZE_I(array)
#    define BOOST_PP_ARRAY_SIZE_I(array) BOOST_PP_ARRAY_SIZE_II array
#    define BOOST_PP_ARRAY_SIZE_II(size, data) size
# endif
#
# endif
