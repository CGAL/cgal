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
# ifndef BOOST_PREPROCESSOR_DETAIL_SPLIT_HPP
# define BOOST_PREPROCESSOR_DETAIL_SPLIT_HPP
#
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_SPLIT */
#
# if BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#    define BOOST_PP_SPLIT(n, im) BOOST_PP_SPLIT_I((n, im))
#    define BOOST_PP_SPLIT_I(par) BOOST_PP_SPLIT_II ## par
#    define BOOST_PP_SPLIT_II(n, a, b) BOOST_PP_SPLIT_ ## n(a, b)
# elif BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MSVC()
#    define BOOST_PP_SPLIT(n, im) BOOST_PP_SPLIT_I(n((im)))
#    define BOOST_PP_SPLIT_I(n) BOOST_PP_SPLIT_ID(BOOST_PP_SPLIT_II_ ## n)
#    define BOOST_PP_SPLIT_II_0(s) BOOST_PP_SPLIT_ID(BOOST_PP_SPLIT_0 s)
#    define BOOST_PP_SPLIT_II_1(s) BOOST_PP_SPLIT_ID(BOOST_PP_SPLIT_1 s)
#    define BOOST_PP_SPLIT_ID(id) id
# else
#    define BOOST_PP_SPLIT(n, im) BOOST_PP_SPLIT_I(n)(im)
#    define BOOST_PP_SPLIT_I(n) BOOST_PP_SPLIT_ ## n
# endif
#
# define BOOST_PP_SPLIT_0(a, b) a
# define BOOST_PP_SPLIT_1(a, b) b
#
# endif
