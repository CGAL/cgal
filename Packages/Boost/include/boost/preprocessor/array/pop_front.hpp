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
# ifndef BOOST_PREPROCESSOR_ARRAY_POP_FRONT_HPP
# define BOOST_PREPROCESSOR_ARRAY_POP_FRONT_HPP
#
# include <boost/preprocessor/arithmetic/dec.hpp>
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <boost/preprocessor/array/elem.hpp>
# include <boost/preprocessor/array/size.hpp>
# include <boost/preprocessor/repetition/enum.hpp>
# include <boost/preprocessor/repetition/deduce_z.hpp>
#
# /* BOOST_PP_ARRAY_POP_FRONT */
#
# define BOOST_PP_ARRAY_POP_FRONT(array) BOOST_PP_ARRAY_POP_FRONT_Z(BOOST_PP_DEDUCE_Z(), array)
#
# /* BOOST_PP_ARRAY_POP_FRONT_Z */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_ARRAY_POP_FRONT_Z(z, array) BOOST_PP_ARRAY_POP_FRONT_I(z, BOOST_PP_ARRAY_SIZE(array), array)
# else
#    define BOOST_PP_ARRAY_POP_FRONT_Z(z, array) BOOST_PP_ARRAY_POP_FRONT_Z_D(z, array)
#    define BOOST_PP_ARRAY_POP_FRONT_Z_D(z, array) BOOST_PP_ARRAY_POP_FRONT_I(z, BOOST_PP_ARRAY_SIZE(array), array)
# endif
#
# define BOOST_PP_ARRAY_POP_FRONT_I(z, size, array) (BOOST_PP_DEC(size), (BOOST_PP_ENUM_ ## z(BOOST_PP_DEC(size), BOOST_PP_ARRAY_POP_FRONT_M, array)))
# define BOOST_PP_ARRAY_POP_FRONT_M(z, n, data) BOOST_PP_ARRAY_ELEM(BOOST_PP_INC(n), data)
#
# endif
