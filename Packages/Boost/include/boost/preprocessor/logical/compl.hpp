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
# ifndef BOOST_PREPROCESSOR_LOGICAL_COMPL_HPP
# define BOOST_PREPROCESSOR_LOGICAL_COMPL_HPP
#
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_COMPL */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#    define BOOST_PP_COMPL(x) BOOST_PP_COMPL_I(x)
# else
#    define BOOST_PP_COMPL(x) BOOST_PP_COMPL_OO((x))
#    define BOOST_PP_COMPL_OO(par) BOOST_PP_COMPL_I ## par
# endif
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MSVC()
#    define BOOST_PP_COMPL_I(x) BOOST_PP_COMPL_ ## x
# else
#    define BOOST_PP_COMPL_I(x) BOOST_PP_COMPL_ID(BOOST_PP_COMPL_ ## x)
#    define BOOST_PP_COMPL_ID(id) id
# endif
#
# define BOOST_PP_COMPL_0 1
# define BOOST_PP_COMPL_1 0
#
# endif
