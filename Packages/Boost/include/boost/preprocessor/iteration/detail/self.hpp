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
# if !defined(BOOST_PP_INDIRECT_SELF)
#    error BOOST_PP_ERROR:  no indirect file to include
# endif
#
# define BOOST_PP_IS_SELFISH 1
#
# include BOOST_PP_INDIRECT_SELF
#
# undef BOOST_PP_IS_SELFISH
# undef BOOST_PP_INDIRECT_SELF
