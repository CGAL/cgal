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
# ifndef BOOST_PREPROCESSOR_REPETITION_ENUM_PARAMS_WITH_DEFAULTS_HPP
# define BOOST_PREPROCESSOR_REPETITION_ENUM_PARAMS_WITH_DEFAULTS_HPP
#
# include <boost/preprocessor/config/config.hpp>
# include <boost/preprocessor/repetition/enum_binary_params.hpp>
#
# /* BOOST_PP_ENUM_PARAMS_WITH_DEFAULTS */
#
# define BOOST_PP_ENUM_PARAMS_WITH_DEFAULTS(count, param, def) BOOST_PP_ENUM_BINARY_PARAMS(count, param, = def)
#
# endif
