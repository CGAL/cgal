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
# ifndef BOOST_PREPROCESSOR_SEQ_REMOVE_HPP
# define BOOST_PREPROCESSOR_SEQ_REMOVE_HPP
#
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/seq/first_n.hpp>
# include <boost/preprocessor/seq/rest_n.hpp>
#
# /* BOOST_PP_SEQ_REMOVE */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_EDG()
#    define BOOST_PP_SEQ_REMOVE(seq, i) BOOST_PP_SEQ_FIRST_N(i, seq) BOOST_PP_SEQ_REST_N(BOOST_PP_INC(i), seq)
# else
#    define BOOST_PP_SEQ_REMOVE(seq, i) BOOST_PP_SEQ_REMOVE_I(seq, i)
#    define BOOST_PP_SEQ_REMOVE_I(seq, i) BOOST_PP_SEQ_FIRST_N(i, seq) BOOST_PP_SEQ_REST_N(BOOST_PP_INC(i), seq)
# endif
#
# endif
