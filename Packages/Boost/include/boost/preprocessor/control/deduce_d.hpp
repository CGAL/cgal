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
# ifndef BOOST_PREPROCESSOR_CONTROL_DEDUCE_D_HPP
# define BOOST_PREPROCESSOR_CONTROL_DEDUCE_D_HPP
#
# include <boost/preprocessor/control/while.hpp>
# include <boost/preprocessor/detail/auto_rec.hpp>
#
# /* BOOST_PP_DEDUCE_D */
#
# define BOOST_PP_DEDUCE_D() BOOST_PP_AUTO_REC(BOOST_PP_WHILE_P, 256)
#
# endif
