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
# ifndef BOOST_PREPROCESSOR_LOGICAL_BITXOR_HPP
# define BOOST_PREPROCESSOR_LOGICAL_BITXOR_HPP
#
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_BITXOR */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#    define BOOST_PP_BITXOR(x, y) BOOST_PP_BITXOR_I(x, y)
# else
#    define BOOST_PP_BITXOR(x, y) BOOST_PP_BITXOR_OO((x, y))
#    define BOOST_PP_BITXOR_OO(par) BOOST_PP_BITXOR_I ## par
# endif
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MSVC()
#    define BOOST_PP_BITXOR_I(x, y) BOOST_PP_BITXOR_ ## x ## y
# else
#    define BOOST_PP_BITXOR_I(x, y) BOOST_PP_BITXOR_ID(BOOST_PP_BITXOR_ ## x ## y)
#    define BOOST_PP_BITXOR_ID(id) id
# endif
#
# define BOOST_PP_BITXOR_00 0
# define BOOST_PP_BITXOR_01 1
# define BOOST_PP_BITXOR_10 1
# define BOOST_PP_BITXOR_11 0
#
# endif
