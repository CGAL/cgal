# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2003.  Permission to copy, use,        *
#  *     modify, sell, and distribute this software is granted provided       *
#  *     this copyright notice appears in all copies.  This software is       *
#  *     provided "as is" without express or implied warranty, and with       *
#  *     no claim at to its suitability for any purpose.                      *
#  *                                                                          *
#  ************************************************************************** */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_FACILITIES_IS_1_HPP
# define BOOST_PREPROCESSOR_FACILITIES_IS_1_HPP
#
# include <boost/preprocessor/cat.hpp>
# include <boost/preprocessor/facilities/is_empty.hpp>
#
# /* BOOST_PP_IS_1 */
#
# define BOOST_PP_IS_1(x) BOOST_PP_IS_EMPTY(BOOST_PP_CAT(BOOST_PP_IS_1_HELPER_, x))
# define BOOST_PP_IS_1_HELPER_1
#
# endif
