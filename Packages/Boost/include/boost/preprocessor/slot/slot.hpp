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
# ifndef BOOST_PREPROCESSOR_SLOT_SLOT_HPP
# define BOOST_PREPROCESSOR_SLOT_SLOT_HPP
#
# include <boost/preprocessor/cat.hpp>
# include <boost/preprocessor/slot/detail/def.hpp>
#
# /* BOOST_PP_ASSIGN_SLOT */
#
# define BOOST_PP_ASSIGN_SLOT(i) BOOST_PP_CAT(BOOST_PP_ASSIGN_SLOT_, i)
#
# define BOOST_PP_ASSIGN_SLOT_1 <boost/preprocessor/slot/detail/slot1.hpp>
# define BOOST_PP_ASSIGN_SLOT_2 <boost/preprocessor/slot/detail/slot2.hpp>
# define BOOST_PP_ASSIGN_SLOT_3 <boost/preprocessor/slot/detail/slot3.hpp>
# define BOOST_PP_ASSIGN_SLOT_4 <boost/preprocessor/slot/detail/slot4.hpp>
# define BOOST_PP_ASSIGN_SLOT_5 <boost/preprocessor/slot/detail/slot5.hpp>
#
# /* BOOST_PP_SLOT */
#
# define BOOST_PP_SLOT(i) BOOST_PP_CAT(BOOST_PP_SLOT_, i)()
#
# endif
