# /* Copyright (C) 2001
#  * Housemarque Oy
#  * http://www.housemarque.com
#  *
#  * Permission to copy, use, modify, sell and distribute this software is
#  * granted provided this copyright notice appears in all copies. This
#  * software is provided "as is" without express or implied warranty, and
#  * with no claim as to its suitability for any purpose.
#  */
#
# /* Revised by Paul Mensonides (2002) */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_STRINGIZE_HPP
# define BOOST_PREPROCESSOR_STRINGIZE_HPP
#
# include <boost/preprocessor/config/config.hpp>
#
# /* BOOST_PP_STRINGIZE */
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_MWCC()
#    define BOOST_PP_STRINGIZE(text) BOOST_PP_STRINGIZE_I(text)
# else
#    define BOOST_PP_STRINGIZE(text) BOOST_PP_STRINGIZE_OO((text))
#    define BOOST_PP_STRINGIZE_OO(par) BOOST_PP_STRINGIZE_I ## par
# endif
#
# define BOOST_PP_STRINGIZE_I(text) #text
#
# endif
