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
# ifndef BOOST_PREPROCESSOR_SEQ_POP_FRONT_HPP
# define BOOST_PREPROCESSOR_SEQ_POP_FRONT_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/seq/seq.hpp>
#
# /* BOOST_PP_SEQ_POP_FRONT */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_POP_FRONT(seq) BOOST_PP_SEQ_TAIL(seq)
# else
#    define BOOST_PP_SEQ_POP_FRONT(seq) BOOST_PP_SEQ_POP_FRONT_I(seq)
#    define BOOST_PP_SEQ_POP_FRONT_I(seq) BOOST_PP_SEQ_TAIL(seq)
# endif
#
# endif
