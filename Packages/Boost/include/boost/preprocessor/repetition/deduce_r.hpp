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
# ifndef BOOST_PREPROCESSOR_REPETITION_DEDUCE_R_HPP
# define BOOST_PREPROCESSOR_REPETITION_DEDUCE_R_HPP
#
# include <boost/preprocessor/detail/auto_rec.hpp>
# include <boost/preprocessor/repetition/for.hpp>
#
# /* BOOST_PP_DEDUCE_R */
#
# define BOOST_PP_DEDUCE_R() BOOST_PP_AUTO_REC(BOOST_PP_FOR_P, 256)
#
# endif
